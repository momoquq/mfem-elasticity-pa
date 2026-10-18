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
namespace internal
{

// Add elasticity action
template <int D1D, int Q1D>
void My3DAddMultPA_(
    const int nDofs,
    const mfem::FiniteElementSpace &fespace,
    const mfem::CoefficientVector &lambda,
    const mfem::CoefficientVector &mu,
    const mfem::GeometricFactors &geom,
    const mfem::DofToQuad &maps,
    const mfem::Vector &x,
    mfem::QuadratureFunction &QVec,
    mfem::Vector &y);

// Assemble element diagonal
template <int D1D, int Q1D>
void My3DAssembleDiagonalPA_(
    const int nDofs,
    const mfem::CoefficientVector &lambda,
    const mfem::CoefficientVector &mu,
    const mfem::GeometricFactors &geom,
    const mfem::DofToQuad &maps,
    mfem::QuadratureFunction &QVec,
    mfem::Vector &diag);


} // namespace internal
} // namespace mfem
