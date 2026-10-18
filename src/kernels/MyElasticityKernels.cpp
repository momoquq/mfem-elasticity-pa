// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
//
// Project modifications Copyright (c) 2026, MFEM Elasticity PA contributors.

#include "MyElasticityKernels.hpp"

#include "mfem/config/config.hpp"
#include "mfem/general/array.hpp"
#include "mfem/general/forall.hpp"
#include "mfem/linalg/dtensor.hpp"
#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/tensor.hpp"
#include "mfem/fem/quadinterpolator.hpp"
#include "mfem/fem/bilininteg.hpp"
#include "mfem/fem/coefficient.hpp"
#include "mfem/fem/qfunction.hpp"

#include "mfem/fem/kernels.hpp"


#include "mfem.hpp"

namespace mfem
{

namespace internal
{


template <int D1D, int Q1D>
void My3DAddMultPA_(
    const int nDofs,
    const FiniteElementSpace &fespace,
    const CoefficientVector &lambda,
    const CoefficientVector &mu,
    const GeometricFactors &geom,
    const DofToQuad &maps,
    const Vector &x,
    QuadratureFunction &QVec,
    Vector &y)
{

    

    static constexpr int d = 3;
    static constexpr int qSize = d;
    
    // Storage padding
    static constexpr int SIMD_WIDTH = 4;
    static constexpr int ALIGNED_Q1D = ((Q1D + SIMD_WIDTH - 1) / SIMD_WIDTH) * SIMD_WIDTH;
    static constexpr int ALIGNED_D1D = ((D1D + SIMD_WIDTH - 1) / SIMD_WIDTH) * SIMD_WIDTH;

    const auto &ir = *maps.IntRule;
    DofToQuad tensorMaps = maps.FE->GetDofToQuad(ir, DofToQuad::TENSOR);
    const int numEls = fespace.GetNE();
    if (numEls == 0) { return; }

    const int numPoints = ir.GetNPoints();

    // Data views
    const auto lamDev = Reshape(lambda.Read(), numPoints, numEls);
    const auto muDev = Reshape(mu.Read(), numPoints, numEls);
    const auto J_glob = Reshape(geom.J.Read(), numPoints, d, d, numEls);
    const real_t *ipWeights = ir.GetWeights().Read();
    auto Gview = Reshape(tensorMaps.G.Read(), Q1D, D1D);
    auto Bview = Reshape(tensorMaps.B.Read(), Q1D, D1D);

    auto xDev = Reshape(x.Read(), D1D, D1D, D1D, d, numEls);
    auto yDev = Reshape(y.ReadWrite(), D1D, D1D, D1D, qSize, numEls);

    mfem::forall_2D(numEls, Q1D, Q1D, [=] MFEM_HOST_DEVICE (int e)
    {
        // Gradient storage.
        alignas(64) MFEM_SHARED real_t grad_vec[d][d][ALIGNED_Q1D][ALIGNED_Q1D][ALIGNED_Q1D];
        
        // Inverse Jacobian storage
        alignas(64) MFEM_SHARED real_t J_inv[3][3];
        alignas(64) MFEM_SHARED real_t detJ_inv;
        
        // Compute the inverse Jacobian
        {
            const real_t J00 = J_glob(0, 0, 0, e), J01 = J_glob(0, 0, 1, e), J02 = J_glob(0, 0, 2, e);
            const real_t J10 = J_glob(0, 1, 0, e), J11 = J_glob(0, 1, 1, e), J12 = J_glob(0, 1, 2, e);
            const real_t J20 = J_glob(0, 2, 0, e), J21 = J_glob(0, 2, 1, e), J22 = J_glob(0, 2, 2, e);
            
            const real_t detJ = J00*(J11*J22 - J12*J21) - J01*(J10*J22 - J12*J20) + J02*(J10*J21 - J11*J20);
            detJ_inv = 1.0 / detJ;
            
            J_inv[0][0] = (J11*J22 - J12*J21) * detJ_inv;
            J_inv[0][1] = (J02*J21 - J01*J22) * detJ_inv;
            J_inv[0][2] = (J01*J12 - J02*J11) * detJ_inv;
            J_inv[1][0] = (J12*J20 - J10*J22) * detJ_inv;
            J_inv[1][1] = (J00*J22 - J02*J20) * detJ_inv;
            J_inv[1][2] = (J02*J10 - J00*J12) * detJ_inv;
            J_inv[2][0] = (J10*J21 - J11*J20) * detJ_inv;
            J_inv[2][1] = (J01*J20 - J00*J21) * detJ_inv;
            J_inv[2][2] = (J00*J11 - J01*J10) * detJ_inv;
        }

        // Load basis matrices
        constexpr int MQ1 = DofQuadLimits::MAX_INTERP_1D;
        constexpr int MD1 = DofQuadLimits::MAX_INTERP_1D;
        alignas(64) MFEM_SHARED real_t BG[2][MQ1*MD1];
        
        DeviceMatrix B(BG[0], D1D, Q1D);
        DeviceMatrix G(BG[1], D1D, Q1D);
        kernels::internal::LoadBG<MD1,MQ1>(D1D, Q1D, Bview, Gview, BG);

        // Forward workspaces
        alignas(64) MFEM_SHARED real_t sm0[2][ALIGNED_Q1D*ALIGNED_D1D];
        alignas(64) MFEM_SHARED real_t sm1[3][ALIGNED_Q1D*ALIGNED_Q1D*ALIGNED_D1D];

        // Forward sum factorization
        for (int c = 0; c < d; ++c)
        {
            // Process a z-slice
            for (int iz = 0; iz < D1D; ++iz)
            {
                // Contract x
                MFEM_FOREACH_THREAD(iy,y,D1D) 
                {
                    #pragma omp simd aligned(sm0:64)
                    for (int qx = 0; qx < Q1D; ++qx) 
                    {
                        real_t u = 0.0, v = 0.0;
                        
                        #pragma unroll(4)
                        for (int ix = 0; ix < D1D; ++ix) {
                            const real_t val = xDev(ix,iy,iz,c,e);
                            u += val * B(ix,qx);
                            v += val * G(ix,qx);
                        }
                        
                        const int idx = iy*ALIGNED_Q1D + qx;
                        sm0[0][idx] = u;
                        sm0[1][idx] = v;
                    }
                }
                MFEM_SYNC_THREAD;
                
                // Contract y
                MFEM_FOREACH_THREAD(qy,y,Q1D) 
                {
                    #pragma omp simd aligned(sm1:64)
                    for (int qx = 0; qx < Q1D; ++qx) 
                    {
                        real_t u = 0.0, v = 0.0, w = 0.0;
                        
                        #pragma unroll(4)
                        for (int iy = 0; iy < D1D; ++iy) {
                            const int idx_in = iy*ALIGNED_Q1D + qx;
                            u += sm0[1][idx_in] * B(iy,qy);  // dxdu * B
                            v += sm0[0][idx_in] * G(iy,qy);  // u * dyB
                            w += sm0[0][idx_in] * B(iy,qy);  // u * B
                        }
                        
                        const int idx_out = iz*ALIGNED_Q1D*ALIGNED_Q1D + qy*ALIGNED_Q1D + qx;
                        sm1[0][idx_out] = u;  // dx
                        sm1[1][idx_out] = v;  // dy
                        sm1[2][idx_out] = w;  // value
                    }
                }
                MFEM_SYNC_THREAD;
            }
            
            // Contract z
            MFEM_FOREACH_THREAD(qz,z,Q1D) 
            {
                for (int qy = 0; qy < Q1D; ++qy) 
                {
                    #pragma omp simd aligned(grad_vec:64)
                    for (int qx = 0; qx < Q1D; ++qx) 
                    {
                        real_t grad_xi[3] = {0.0, 0.0, 0.0};
                        
                        #pragma unroll(4)
                        for (int iz = 0; iz < D1D; ++iz) {
                            const int idx = iz*ALIGNED_Q1D*ALIGNED_Q1D + qy*ALIGNED_Q1D + qx;
                            grad_xi[0] += sm1[0][idx] * B(iz,qz);  // dx
                            grad_xi[1] += sm1[1][idx] * B(iz,qz);  // dy
                            grad_xi[2] += sm1[2][idx] * G(iz,qz);  // dz
                        }
                        
                        // Transform gradients
                        #pragma unroll(3)
                        for (int i = 0; i < 3; ++i) {
                            grad_vec[c][i][qz][qy][qx] = J_inv[0][i]*grad_xi[0] +
                                                        J_inv[1][i]*grad_xi[1] +
                                                        J_inv[2][i]*grad_xi[2];
                        }
                    }
                }
            }
            MFEM_SYNC_THREAD;
        }

        // Stress storage
        alignas(64) MFEM_SHARED real_t stress[6][ALIGNED_Q1D][ALIGNED_Q1D][ALIGNED_Q1D];
        
        MFEM_FOREACH_THREAD(qz,z,Q1D) 
        {
            for (int qy = 0; qy < Q1D; ++qy) 
            {
                #pragma omp simd aligned(stress,grad_vec:64)
                for (int qx = 0; qx < Q1D; ++qx) 
                {
                    const int p = qx + Q1D * (qy + Q1D * qz);
                    
                    // Divergence
                    const real_t div = grad_vec[0][0][qz][qy][qx] + 
                                        grad_vec[1][1][qz][qy][qx] + 
                                        grad_vec[2][2][qz][qy][qx];
                    
                    const real_t w = ipWeights[p] / detJ_inv;
                    const real_t lam_w = lamDev(p, e) * w;
                    const real_t mu_w = muDev(p, e) * w;
                    
                    // Compute weighted stress
                    stress[0][qz][qy][qx] = lam_w * div + 2.0 * mu_w * grad_vec[0][0][qz][qy][qx];  // sigma_00
                    stress[1][qz][qy][qx] = lam_w * div + 2.0 * mu_w * grad_vec[1][1][qz][qy][qx];  // sigma_11
                    stress[2][qz][qy][qx] = lam_w * div + 2.0 * mu_w * grad_vec[2][2][qz][qy][qx];  // sigma_22
                    stress[3][qz][qy][qx] = mu_w * (grad_vec[0][1][qz][qy][qx] + grad_vec[1][0][qz][qy][qx]);  // sigma_01 = sigma_10
                    stress[4][qz][qy][qx] = mu_w * (grad_vec[0][2][qz][qy][qx] + grad_vec[2][0][qz][qy][qx]);  // sigma_02 = sigma_20
                    stress[5][qz][qy][qx] = mu_w * (grad_vec[1][2][qz][qy][qx] + grad_vec[2][1][qz][qy][qx]);  // sigma_12 = sigma_21
                }
            }
        }
        MFEM_SYNC_THREAD;

        // Backward sum factorization
        for (int q = 0; q < qSize; ++q) 
        {
            // Displacement component q
            alignas(64) MFEM_SHARED real_t Q_vec_q[3][ALIGNED_Q1D][ALIGNED_Q1D][ALIGNED_Q1D];
            
            MFEM_FOREACH_THREAD(qz,z,Q1D) 
            {
                for (int qy = 0; qy < Q1D; ++qy) 
                {
                    #pragma omp simd aligned(Q_vec_q,stress:64)
                    for (int qx = 0; qx < Q1D; ++qx) 
                    {
                        // Reference flux
                        #pragma unroll(3)
                        for (int m = 0; m < 3; ++m) {
                            real_t sum = 0.0;
                            
                            if (q == 0) {
                                // Stress row 0
                                sum = stress[0][qz][qy][qx] * J_inv[m][0] +   // sigma_00
                                        stress[3][qz][qy][qx] * J_inv[m][1] +   // sigma_01
                                        stress[4][qz][qy][qx] * J_inv[m][2];    // sigma_02
                            } else if (q == 1) {
                                // Stress row 1
                                sum = stress[3][qz][qy][qx] * J_inv[m][0] +   // sigma_10 = sigma_01
                                        stress[1][qz][qy][qx] * J_inv[m][1] +   // sigma_11
                                        stress[5][qz][qy][qx] * J_inv[m][2];    // sigma_12
                            } else { // q == 2
                                // Stress row 2
                                sum = stress[4][qz][qy][qx] * J_inv[m][0] +   // sigma_20 = sigma_02
                                        stress[5][qz][qy][qx] * J_inv[m][1] +   // sigma_21 = sigma_12
                                        stress[2][qz][qy][qx] * J_inv[m][2];    // sigma_22
                            }
                            
                            Q_vec_q[m][qz][qy][qx] = sum;
                        }
                    }
                }
            }
            MFEM_SYNC_THREAD;
            
            // Process a z-slice
            for (int iz = 0; iz < D1D; ++iz) 
            {
                alignas(64) real_t tmpZ[3][ALIGNED_Q1D][ALIGNED_Q1D];
                
                // Contract z
                MFEM_FOREACH_THREAD(qy,y,Q1D) 
                {
                    for (int qx = 0; qx < Q1D; ++qx) 
                    {
                        #pragma unroll(3)
                        for (int m = 0; m < 3; ++m) {
                            real_t sum = 0.0;
                            const real_t* g_ptr = (m == 2) ? &Gview(0, iz) : &Bview(0, iz);
                            
                            #pragma unroll(8)
                            for (int qz = 0; qz < Q1D; ++qz) {
                                sum += Q_vec_q[m][qz][qy][qx] * g_ptr[qz];
                            }
                            tmpZ[m][qy][qx] = sum;
                        }
                    }
                }
                MFEM_SYNC_THREAD;
                
                // Contract y
                alignas(64) real_t tmpY[3][ALIGNED_D1D][ALIGNED_Q1D];
                
                MFEM_FOREACH_THREAD(iy,y,D1D) 
                {
                    #pragma omp simd aligned(tmpY,tmpZ:64)
                    for (int qx = 0; qx < Q1D; ++qx) 
                    {
                        #pragma unroll(3)
                        for (int m = 0; m < 3; ++m) {
                            real_t sum = 0.0;
                            const real_t* g_ptr = (m == 1) ? &Gview(0, iy) : &Bview(0, iy);
                            
                            #pragma unroll(8)
                            for (int qy = 0; qy < Q1D; ++qy) {
                                sum += tmpZ[m][qy][qx] * g_ptr[qy];
                            }
                            tmpY[m][iy][qx] = sum;
                        }
                    }
                }
                MFEM_SYNC_THREAD;
                
                // Contract x
                MFEM_FOREACH_THREAD(iy,y,D1D) 
                {
                    for (int ix = 0; ix < D1D; ++ix) 
                    {
                        real_t result = 0.0;
                        
                        #pragma unroll(3)
                        for (int m = 0; m < 3; ++m) {
                            const real_t* g_ptr = (m == 0) ? &Gview(0, ix) : &Bview(0, ix);
                            real_t sum = 0.0;
                            
                            #pragma omp simd reduction(+:sum)
                            for (int qx = 0; qx < Q1D; ++qx) {
                                sum += tmpY[m][iy][qx] * g_ptr[qx];
                            }
                            result += sum;
                        }
                        
                        yDev(ix, iy, iz, q, e) += result;
                    }
                }
                MFEM_SYNC_THREAD;
            }
        }
    });


}

template <int D1D, int Q1D>
void My3DAssembleDiagonalPA_(
    const int nDofs,
    const mfem::CoefficientVector &lambda,
    const mfem::CoefficientVector &mu,
    const mfem::GeometricFactors &geom,
    const mfem::DofToQuad &maps,
    mfem::QuadratureFunction &QVec,
    mfem::Vector &diag)
{


    // Element quadrature rule
    const auto &ir = QVec.GetIntRule(0);
    static constexpr int d = 3;
    const int numPoints = ir.GetNPoints();
    const int numEls = lambda.Size()/numPoints;
    const auto lamDev = Reshape(lambda.Read(), numPoints, numEls);
    const auto muDev = Reshape(mu.Read(), numPoints, numEls);
    const auto J = Reshape(geom.J.Read(), numPoints, d, d, numEls);
    auto Q = Reshape(QVec.ReadWrite(), numPoints, d,d, d, numEls);
    const real_t *ipWeights = ir.GetWeights().Read();
    mfem::forall_2D(numEls, numPoints,1, [=] MFEM_HOST_DEVICE (int e)
    {
        MFEM_FOREACH_THREAD(p, x,numPoints)
        {
            auto invJ = inv(make_tensor<d, d>(
                [&](int i, int j) { return J(p, i, j, e); }));
            const real_t w = ipWeights[p] /det(invJ);
            for (int n = 0; n < d; n++)
            {
                for (int m = 0; m < d; m++)
                {
                    for (int q = 0; q < d; q++)
                    {
                        // Symmetric gradient contraction
                        real_t contraction = 0.;
                        for (int a = 0; a < d; a++)
                        {
                            for (int b = 0; b < d; b++)
                            {
                                contraction += ((a == q)*invJ(m,b) + (b==q)*invJ(m,a))*((a == q)
                                    *invJ(n, b) + (b==q)*invJ(n,a));
                            }
                        }
                        // Weighted elasticity tensor
                        Q(p,m,n,q,e) = w*(lamDev(p, e)*invJ(m,q)*invJ(n,q)
                                    + 0.5*muDev(p, e)*contraction);
                    }
                }
            }
        }
    });

    // Reduce to the E-vector
    DofToQuad tensorMaps = maps.FE->GetDofToQuad(*maps.IntRule, DofToQuad::TENSOR);

    //  const int D1D = tensorMaps.ndof;
    //  const int Q1D = tensorMaps.nqpt;

    auto G_1d = Reshape(tensorMaps.G.Read(), Q1D, D1D);
    auto B_1d = Reshape(tensorMaps.B.Read(), Q1D, D1D);
    const auto QRead = Reshape(QVec.Read(), Q1D, Q1D, Q1D, d, d, d, numEls);
    auto diagDev = Reshape(diag.Write(), D1D, D1D, D1D, d, numEls);

    mfem::forall_2D(numEls, d, nDofs, [=] MFEM_HOST_DEVICE (int e)
    {
        MFEM_FOREACH_THREAD(q, x, d)
        {
            for (int n = 0; n < 3; ++n) {
                for (int m = 0; m < 3; ++m) {

                    // Diagonal workspaces
                    real_t T1[D1D][Q1D][Q1D];
                    real_t T2[D1D][D1D][Q1D];

                    // Contract x
                    for (int dx = 0; dx < D1D; ++dx) {
                        for (int qz = 0; qz < Q1D; ++qz) {
                            for (int qy = 0; qy < Q1D; ++qy) {
                                real_t sum = 0.0;
                                for (int qx = 0; qx < Q1D; ++qx) {
                                    const double Bx = B_1d(qx, dx);
                                    const double Gx = G_1d(qx, dx);
                                    const double op_m = (m == 0) ? Gx : Bx;
                                    const double op_n = (n == 0) ? Gx : Bx;
                                    sum += QRead(qx, qy, qz, m, n, q, e) * op_m * op_n;
                                }
                                T1[dx][qz][qy] = sum;
                            }
                        }
                    }

                    // Contract y
                    for (int dy = 0; dy < D1D; ++dy) {
                        for (int dx = 0; dx < D1D; ++dx) {
                            for (int qz = 0; qz < Q1D; ++qz) {
                                real_t sum = 0.0;
                                for (int qy = 0; qy < Q1D; ++qy) {
                                    const double By = B_1d(qy, dy);
                                    const double Gy = G_1d(qy, dy);
                                    const double op_m = (m == 1) ? Gy : By;
                                    const double op_n = (n == 1) ? Gy : By;
                                    sum += T1[dx][qz][qy] * op_m * op_n;
                                }
                                T2[dy][dx][qz] = sum;
                            }
                        }
                    }

                    // Contract z
                    for (int dz = 0; dz < D1D; ++dz) {
                        for (int dy = 0; dy < D1D; ++dy) {
                            for (int dx = 0; dx < D1D; ++dx) {
                                real_t final_sum = 0.0;
                                for (int qz = 0; qz < Q1D; ++qz) {
                                    const double Bz = B_1d(qz, dz);
                                    const double Gz = G_1d(qz, dz);
                                    const double op_m = (m == 2) ? Gz : Bz;
                                    const double op_n = (n == 2) ? Gz : Bz;
                                    final_sum += T2[dy][dx][qz] * op_m * op_n;
                                }
                                diagDev(dx, dy, dz, q, e) += final_sum;
                            }
                        }
                    }
                }
            }
        }
    });


}

// Explicit action kernel instantiations.

// p=1: <2, 3>
template void My3DAddMultPA_<2, 3>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=2: <3, 4>
template void My3DAddMultPA_<3, 4>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=3: <4, 5>
template void My3DAddMultPA_<4, 5>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=4: <5, 6>
template void My3DAddMultPA_<5, 6>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=5: <6, 7>
template void My3DAddMultPA_<6, 7>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=6: <7, 8>
template void My3DAddMultPA_<7, 8>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=7: <8, 9>
template void My3DAddMultPA_<8, 9>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=8: <9, 10>
template void My3DAddMultPA_<9, 10>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=9: <10, 11>
template void My3DAddMultPA_<10, 11>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=10: <11, 12>
template void My3DAddMultPA_<11, 12>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=11: <12, 13>
template void My3DAddMultPA_<12, 13>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=12: <13, 14>
template void My3DAddMultPA_<13, 14>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=13: <14, 15>
template void My3DAddMultPA_<14, 15>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=14: <15, 16>
template void My3DAddMultPA_<15, 16>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=15: <16, 17>
template void My3DAddMultPA_<16, 17>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);

// p=16: <17, 18>
template void My3DAddMultPA_<17, 18>(
    const int nDofs, const FiniteElementSpace &fespace,
    const CoefficientVector &lambda, const CoefficientVector &mu,
    const GeometricFactors &geom, const DofToQuad &maps,
    const Vector &x, QuadratureFunction &QVec, Vector &y);


// Explicit diagonal kernel instantiations.

// p=1: <2, 3>
template void My3DAssembleDiagonalPA_<2, 3>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=2: <3, 4>
template void My3DAssembleDiagonalPA_<3, 4>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=3: <4, 5>
template void My3DAssembleDiagonalPA_<4, 5>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=4: <5, 6>
template void My3DAssembleDiagonalPA_<5, 6>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=5: <6, 7>
template void My3DAssembleDiagonalPA_<6, 7>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=6: <7, 8>
template void My3DAssembleDiagonalPA_<7, 8>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=7: <8, 9>
template void My3DAssembleDiagonalPA_<8, 9>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=8: <9, 10>
template void My3DAssembleDiagonalPA_<9, 10>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=9: <10, 11>
template void My3DAssembleDiagonalPA_<10, 11>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=10: <11, 12>
template void My3DAssembleDiagonalPA_<11, 12>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=11: <12, 13>
template void My3DAssembleDiagonalPA_<12, 13>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=12: <13, 14>
template void My3DAssembleDiagonalPA_<13, 14>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);

// p=13: <14, 15>
template void My3DAssembleDiagonalPA_<14, 15>(
    const int nDofs, const CoefficientVector &lambda,
    const CoefficientVector &mu, const GeometricFactors &geom,
    const DofToQuad &maps, QuadratureFunction &QVec, Vector &diag);


} // namespace internal

} // namespace mfem
