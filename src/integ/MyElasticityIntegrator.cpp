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

#include "MyElasticityIntegrator.hpp"
#include "../kernels/MyElasticityKernels.hpp"

namespace mfem
{

   void MyElasticityIntegrator::AssembleDiagonalPA(Vector &diag)
   {

      q_vec->SetVDim(vdim*vdim*vdim*vdim);

      // Check dimension
      MFEM_ASSERT(vdim == 3, "This function only supports 3D problems.");
   
      // Tensor basis map
      const mfem::DofToQuad &tensorMaps = maps->FE->GetDofToQuad(*maps->IntRule, mfem::DofToQuad::TENSOR);
   
      // Basis and quadrature sizes
      const int D1D = tensorMaps.ndof;
      const int Q1D = tensorMaps.nqpt;

      // Encode kernel ID
      const int id = (D1D << 4) | Q1D;

      // Kernel dispatch
      switch (id)
      {
         // p=1: D1D=2, Q1D=3
         case 0x23:
         {
            internal::My3DAssembleDiagonalPA_<2, 3>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=2: D1D=3, Q1D=4
         case 0x34:
         {
            internal::My3DAssembleDiagonalPA_<3, 4>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=3: D1D=4, Q1D=5
         case 0x45:
         {
            internal::My3DAssembleDiagonalPA_<4, 5>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=4: D1D=5, Q1D=6
         case 0x56:
         {
            internal::My3DAssembleDiagonalPA_<5, 6>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=5: D1D=6, Q1D=7
         case 0x67:
         {
            internal::My3DAssembleDiagonalPA_<6, 7>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=6: D1D=7, Q1D=8
         case 0x78:
         {
            internal::My3DAssembleDiagonalPA_<7, 8>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=7: D1D=8, Q1D=9
         case 0x89:
         {
            internal::My3DAssembleDiagonalPA_<8, 9>(ndofs, *lambda_quad, *mu_quad,
                                                     *geom, *maps, *q_vec, diag);
            break;
         }
         // p=8: D1D=9, Q1D=10
         case 0x9A:
         {
            internal::My3DAssembleDiagonalPA_<9, 10>(ndofs, *lambda_quad, *mu_quad,
                                                      *geom, *maps, *q_vec, diag);
            break;
         }
         // p=9: D1D=10, Q1D=11
         case 0xAB:
         {
            internal::My3DAssembleDiagonalPA_<10, 11>(ndofs, *lambda_quad, *mu_quad,
                                                       *geom, *maps, *q_vec, diag);
            break;
         }
         // p=10: D1D=11, Q1D=12
         case 0xBC:
         {
            internal::My3DAssembleDiagonalPA_<11, 12>(ndofs, *lambda_quad, *mu_quad,
                                                       *geom, *maps, *q_vec, diag);
            break;
         }
         // p=11: D1D=12, Q1D=13
         case 0xCD:
         {
            internal::My3DAssembleDiagonalPA_<12, 13>(ndofs, *lambda_quad, *mu_quad,
                                                       *geom, *maps, *q_vec, diag);
            break;
         }
         // p=12: D1D=13, Q1D=14
         case 0xDE:
         {
            internal::My3DAssembleDiagonalPA_<13, 14>(ndofs, *lambda_quad, *mu_quad,
                                                       *geom, *maps, *q_vec, diag);
            break;
         }
         // p=13: D1D=14, Q1D=15
         case 0xEF:
         {
            internal::My3DAssembleDiagonalPA_<14, 15>(ndofs, *lambda_quad, *mu_quad,
                                                       *geom, *maps, *q_vec, diag);
            break;
         }
         // Unsupported kernel ID
         default:
         {
            // Kernel diagnostics
            std::ostringstream oss;
            oss << "Kernel not implemented for D1D=" << D1D << ", Q1D=" << Q1D
                << " (ID: 0x" << std::hex << id << std::dec << ")";
            MFEM_ABORT(oss.str());
         }
      }
   }

void MyElasticityIntegrator::AddMultPA(const Vector &x, Vector &y) const
{
   // Check dimension
   MFEM_ASSERT(vdim == 3, "This function only supports 3D problems.");

   // Tensor basis map
   const mfem::DofToQuad &tensorMaps = maps->FE->GetDofToQuad(*maps->IntRule, mfem::DofToQuad::TENSOR);

   // Basis and quadrature sizes
   const int D1D = tensorMaps.ndof;
   const int Q1D = tensorMaps.nqpt;

   // Encode kernel ID
   const int id = (D1D << 4) | Q1D;

   // Kernel dispatch
   switch (id)
   {
      // p=1: D1D=2, Q1D=3
      case 0x23:
      {
         internal::My3DAddMultPA_<2, 3>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=2: D1D=3, Q1D=4
      case 0x34:
      {
         internal::My3DAddMultPA_<3, 4>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=3: D1D=4, Q1D=5
      case 0x45:
      {
         internal::My3DAddMultPA_<4, 5>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=4: D1D=5, Q1D=6
      case 0x56:
      {
         internal::My3DAddMultPA_<5, 6>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=5: D1D=6, Q1D=7
      case 0x67:
      {
         internal::My3DAddMultPA_<6, 7>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=6: D1D=7, Q1D=8
      case 0x78:
      {
         internal::My3DAddMultPA_<7, 8>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=7: D1D=8, Q1D=9
      case 0x89:
      {
         internal::My3DAddMultPA_<8, 9>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                        *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=8: D1D=9, Q1D=10
      case 0x9A:
      {
         internal::My3DAddMultPA_<9, 10>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                         *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=9: D1D=10, Q1D=11
      case 0xAB:
      {
         internal::My3DAddMultPA_<10, 11>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                          *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=10: D1D=11, Q1D=12
      case 0xBC:
      {
         internal::My3DAddMultPA_<11, 12>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                          *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=11: D1D=12, Q1D=13
      case 0xCD:
      {
         internal::My3DAddMultPA_<12, 13>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                          *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=12: D1D=13, Q1D=14
      case 0xDE:
      {
         internal::My3DAddMultPA_<13, 14>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                          *geom, *maps, x, *q_vec, y);
         break;
      }
      // p=13: D1D=14, Q1D=15
      case 0xEF:
      {
         internal::My3DAddMultPA_<14, 15>(ndofs, *fespace, *lambda_quad, *mu_quad,
                                          *geom, *maps, x, *q_vec, y);
         break;
      }
      // Unsupported kernel ID
      default:
      {
         std::ostringstream oss;
         oss << "Kernel not implemented for D1D=" << D1D << ", Q1D=" << Q1D
             << " (ID: 0x" << std::hex << id << std::dec << ")";
         MFEM_ABORT(oss.str());
      }
   }
}

void MyElasticityIntegrator::AddMultTransposePA(const Vector &x, Vector &y) const
{
   MyElasticityIntegrator::AddMultPA(x, y); // Symmetric operator
}


// Assemble element matrix
void MyElasticityIntegrator::AssembleElementMatrix(
   const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
{
   int dof = el.GetDof();
   int dim = el.GetDim();
   real_t w, L, M;

   MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

#ifdef MFEM_THREAD_SAFE
   DenseMatrix dshape(dof, dim), gshape(dof, dim), pelmat(dof);
   Vector divshape(dim*dof);
#else
   dshape.SetSize(dof, dim);
   gshape.SetSize(dof, dim);
   pelmat.SetSize(dof);
   divshape.SetSize(dim*dof);
#endif

   elmat.SetSize(dof * dim);

   const IntegrationRule *ir = GetIntegrationRule(el, Trans);
   if (ir == NULL)
   {
      int order = 2 * Trans.OrderGrad(&el);
      ir = &IntRules.Get(el.GetGeomType(), order);
   }

   elmat = 0.0;

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      el.CalcDShape(ip, dshape);

      Trans.SetIntPoint(&ip);
      w = ip.weight * Trans.Weight();
      Mult(dshape, Trans.InverseJacobian(), gshape);
      MultAAt(gshape, pelmat);
      gshape.GradToDiv (divshape);

      M = mu->Eval(Trans, ip);
      if (lambda)
      {
         L = lambda->Eval(Trans, ip);
      }
      else
      {
         L = q_lambda * M;
         M = q_mu * M;
      }

      if (L != 0.0)
      {
         AddMult_a_VVt(L * w, divshape, elmat);
      }

      if (M != 0.0)
      {
         for (int d = 0; d < dim; d++)
         {
            for (int k = 0; k < dof; k++)
               for (int l = 0; l < dof; l++)
               {
                  elmat (dof*d+k, dof*d+l) += (M * w) * pelmat(k, l);
               }
         }
         for (int ii = 0; ii < dim; ii++)
            for (int jj = 0; jj < dim; jj++)
            {
               for (int kk = 0; kk < dof; kk++)
                  for (int ll = 0; ll < dof; ll++)
                  {
                     elmat(dof*ii+kk, dof*jj+ll) +=
                        (M * w) * gshape(kk, jj) * gshape(ll, ii);
                  }
            }
      }
   }
}

// Compute stress
void MyElasticityIntegrator::ComputeElementFlux(
   const mfem::FiniteElement &el, ElementTransformation &Trans,
   Vector &u, const mfem::FiniteElement &fluxelem, Vector &flux,
   bool with_coef, const IntegrationRule *ir)
{
   const int dof = el.GetDof();
   const int dim = el.GetDim();
   const int tdim = dim*(dim+1)/2; // Symmetric tensor size
   real_t L, M;

   MFEM_ASSERT(dim == 2 || dim == 3,
               "dimension is not supported: dim = " << dim);
   MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");
   MFEM_ASSERT(fluxelem.GetMapType() == FiniteElement::VALUE, "");
   MFEM_ASSERT(dynamic_cast<const NodalFiniteElement*>(&fluxelem), "");

#ifdef MFEM_THREAD_SAFE
   DenseMatrix dshape(dof, dim);
#else
   dshape.SetSize(dof, dim);
#endif

   real_t gh_data[9], grad_data[9];
   DenseMatrix gh(gh_data, dim, dim);
   DenseMatrix grad(grad_data, dim, dim);

   if (!ir)
   {
      ir = &fluxelem.GetNodes();
   }
   const int fnd = ir->GetNPoints();
   flux.SetSize(fnd * tdim);

   DenseMatrix loc_data_mat(u.GetData(), dof, dim);
   for (int i = 0; i < fnd; i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      el.CalcDShape(ip, dshape);
      MultAtB(loc_data_mat, dshape, gh);

      Trans.SetIntPoint(&ip);
      Mult(gh, Trans.InverseJacobian(), grad);

      M = mu->Eval(Trans, ip);
      if (lambda)
      {
         L = lambda->Eval(Trans, ip);
      }
      else
      {
         L = q_lambda * M;
         M = q_mu * M;
      }

      // Compute stress
      const real_t M2 = 2.0*M;
      if (dim == 2)
      {
         L *= (grad(0,0) + grad(1,1));
         // xx, yy, xy
         flux(i+fnd*0) = M2*grad(0,0) + L;
         flux(i+fnd*1) = M2*grad(1,1) + L;
         flux(i+fnd*2) = M*(grad(0,1) + grad(1,0));
      }
      else if (dim == 3)
      {
         L *= (grad(0,0) + grad(1,1) + grad(2,2));
         // xx, yy, zz, xy, xz, yz
         flux(i+fnd*0) = M2*grad(0,0) + L;
         flux(i+fnd*1) = M2*grad(1,1) + L;
         flux(i+fnd*2) = M2*grad(2,2) + L;
         flux(i+fnd*3) = M*(grad(0,1) + grad(1,0));
         flux(i+fnd*4) = M*(grad(0,2) + grad(2,0));
         flux(i+fnd*5) = M*(grad(1,2) + grad(2,1));
      }
   }
}

real_t MyElasticityIntegrator::ComputeFluxEnergy(const FiniteElement &fluxelem,
                                               ElementTransformation &Trans,
                                               Vector &flux, Vector *d_energy)
{
   const int dof = fluxelem.GetDof();
   const int dim = fluxelem.GetDim();
   const int tdim = dim*(dim+1)/2; // Symmetric tensor size
   real_t L, M;

   // Check flux inputs
   MFEM_ASSERT(d_energy == NULL, "anisotropic estimates are not supported");
   MFEM_ASSERT(flux.Size() == dof*tdim, "invalid 'flux' vector");

#ifndef MFEM_THREAD_SAFE
   shape.SetSize(dof);
#else
   Vector shape(dof);
#endif
   real_t pointstress_data[6];
   Vector pointstress(pointstress_data, tdim);

   // Stress view
   DenseMatrix flux_mat(flux.GetData(), dof, tdim);

   // Select quadrature rule

   const IntegrationRule *ir = GetIntegrationRule(fluxelem, Trans);
   if (ir == NULL)
   {
      int order = 2 * Trans.OrderGrad(&fluxelem);
      ir = &IntRules.Get(fluxelem.GetGeomType(), order);
   }

   real_t energy = 0.0;

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Trans.SetIntPoint(&ip);
      fluxelem.CalcPhysShape(Trans, shape);

      flux_mat.MultTranspose(shape, pointstress);

      real_t w = Trans.Weight() * ip.weight;

      M = mu->Eval(Trans, ip);
      if (lambda)
      {
         L = lambda->Eval(Trans, ip);
      }
      else
      {
         L = q_lambda * M;
         M = q_mu * M;
      }

      // Compute strain energy

      real_t pt_e; // Strain energy density
      const real_t *s = pointstress_data;
      if (dim == 2)
      {
         // xx, yy, xy
         const real_t tr_e = (s[0] + s[1])/(2*(M + L));
         L *= tr_e;
         pt_e = (0.25/M)*(s[0]*(s[0] - L) + s[1]*(s[1] - L) + 2*s[2]*s[2]);
      }
      else // (dim == 3)
      {
         // xx, yy, zz, xy, xz, yz
         const real_t tr_e = (s[0] + s[1] + s[2])/(2*M + 3*L);
         L *= tr_e;
         pt_e = (0.25/M)*(s[0]*(s[0] - L) + s[1]*(s[1] - L) + s[2]*(s[2] - L) +
                          2*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
      }

      energy += w * pt_e;
   }
   return energy;
}

// Project material coefficients
void MyElasticityIntegrator::SetUpQuadratureSpaceAndCoefficients(
   const FiniteElementSpace &fes)
{
   if (IntRule == nullptr)
   {
      // Select quadrature rule
      const auto &T = *fes.GetMesh()->GetTypicalElementTransformation();
      int quad_order = 2 * T.OrderGrad(fes.GetTypicalFE());
      IntRule = &IntRules.Get(T.GetGeometryType(), quad_order);
   }

   Mesh &mesh = *fespace->GetMesh();

   q_space.reset(new QuadratureSpace(mesh, *IntRule));
   lambda_quad.reset(new CoefficientVector(lambda, *q_space,
                                           CoefficientStorage::FULL));
   mu_quad.reset(new CoefficientVector(mu, *q_space, CoefficientStorage::FULL));
   q_vec.reset(new QuadratureFunction(*q_space, vdim*vdim));
}

void MyElasticityIntegrator::AssemblePA(const FiniteElementSpace &fes)
{
   MFEM_VERIFY(fes.GetOrdering() == Ordering::byNODES,
               "Elasticity PA only implemented for byNODES ordering.");

   fespace = &fes;
   Mesh &mesh = *fespace->GetMesh();
   MFEM_VERIFY(fespace->GetVDim() == mesh.Dimension(), "");
   vdim = fespace->GetVDim();
   ndofs = fespace->GetTypicalFE()->GetDof();

   SetUpQuadratureSpaceAndCoefficients(fes);

   auto ordering = GetEVectorOrdering(*fespace);
   auto mode = ordering == ElementDofOrdering::NATIVE ? DofToQuad::FULL :
               DofToQuad::LEXICOGRAPHIC_FULL;
   maps = &fespace->GetTypicalFE()->GetDofToQuad(*IntRule, mode);
   geom = mesh.GetGeometricFactors(*IntRule, GeometricFactors::JACOBIANS);
}

}
