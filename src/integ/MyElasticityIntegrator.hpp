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

#pragma once
#include "mfem.hpp"


namespace mfem
{

/// Isotropic elasticity integrator
class MyElasticityIntegrator : public BilinearFormIntegrator
{
protected:
    real_t q_lambda, q_mu;
    Coefficient *lambda, *mu;
 
 private:
 #ifndef MFEM_THREAD_SAFE
    Vector shape;
    DenseMatrix dshape, gshape, pelmat;
    Vector divshape;
 #endif
 
    // PA extension
 
    const DofToQuad *maps;         ///< Borrowed
    const GeometricFactors *geom;  ///< Borrowed
    int vdim, ndofs;
    const FiniteElementSpace *fespace;   ///< Borrowed
 
    std::unique_ptr<QuadratureSpace> q_space;
    /// Quadrature coefficients
    std::unique_ptr<CoefficientVector> lambda_quad, mu_quad;
    /// Workspace vector
    std::unique_ptr<QuadratureFunction> q_vec;
 
    /// Project material coefficients
    void SetUpQuadratureSpaceAndCoefficients(const FiniteElementSpace &fes);
 
 public:
    MyElasticityIntegrator(Coefficient &l, Coefficient &m)
    { lambda = &l; mu = &m; }
    /// Scaled material coefficients
    MyElasticityIntegrator(Coefficient &m, real_t q_l, real_t q_m)
    { lambda = NULL; mu = &m; q_lambda = q_l; q_mu = q_m; }
 
    void AssembleElementMatrix(const FiniteElement &el,
                               ElementTransformation &Tr,
                               DenseMatrix &elmat) override;
 
    using BilinearFormIntegrator::AssemblePA;
    void AssemblePA(const FiniteElementSpace &fes) override;
 
    void AssembleDiagonalPA(Vector &diag) override;
 
    void AddMultPA(const Vector &x, Vector &y) const override;
 
    void AddMultTransposePA(const Vector &x, Vector &y) const override;
 
    /// Compute stress
    void ComputeElementFlux(const FiniteElement &el,
                            ElementTransformation &Trans,
                            Vector &u,
                            const FiniteElement &fluxelem,
                            Vector &flux, bool with_coef = true,
                            const IntegrationRule *ir = NULL) override;
 
    /// Compute strain energy
    real_t ComputeFluxEnergy(const FiniteElement &fluxelem,
                             ElementTransformation &Trans,
                             Vector &flux, Vector *d_energy = NULL) override;

};

} // namespace mfem
