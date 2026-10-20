// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See
// LICENSE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
//
// Project modifications Copyright (c) 2026, MFEM Elasticity PA contributors.

// Elasticity with GMG-PCG

#include "mfem.hpp"
#include "src/integ/MyElasticityIntegrator.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Affine geometry
void checkAffineMesh(Mesh &mesh)
{
   const real_t rtol = 64 * std::numeric_limits<real_t>::epsilon();
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      ElementTransformation *T = mesh.GetElementTransformation(e);
      const IntegrationRule &ir = IntRules.Get(T->GetGeometryType(),
                                               2 * T->OrderJ());
      T->SetIntPoint(&ir.IntPoint(0));
      const DenseMatrix J0(T->Jacobian());
      const int geometry_dofs = mesh.GetNodes()
                                ? mesh.GetNodalFESpace()->GetFE(e)->GetDof()
                                : mesh.GetElement(e)->GetNVertices();
      const real_t tol = rtol * geometry_dofs * J0.MaxMaxNorm();
      for (int q = 0; q < ir.GetNPoints(); ++q)
      {
         T->SetIntPoint(&ir.IntPoint(q));
         const DenseMatrix &J = T->Jacobian();
         for (int i = 0; i < J.Height(); ++i)
         {
            for (int j = 0; j < J.Width(); ++j)
            {
               MFEM_VERIFY(std::abs(J(i, j) - J0(i, j)) <= tol,
                           "Non-affine element " << e
                           << ": optimized PA requires affine geometry.");
            }
         }
      }
   }
}

// Geometric multigrid
class ElasticityMultigrid : public GeometricMultigrid
{
private:
   PWConstCoefficient lambda_func;
   PWConstCoefficient mu_func;
   bool pa;
   bool use_optimized_kernel;
   HypreBoomerAMG* amg;
   bool ownOperator;

public:
   ElasticityMultigrid(ParFiniteElementSpaceHierarchy& fespaces,
                       Array<int>& ess_bdr,
                       const PWConstCoefficient& lambda,
                       const PWConstCoefficient& mu,
                       bool partial_assembly,
                       bool op)
      : GeometricMultigrid(fespaces, ess_bdr),
        lambda_func(lambda),
        mu_func(mu),
        pa(partial_assembly),
        use_optimized_kernel(op)
   {
      ConstructCoarseOperatorAndSolver(fespaces.GetFESpaceAtLevel(0));

      for (int level = 1; level < fespaces.GetNumLevels(); ++level)
      {
         ConstructOperatorAndSmoother(fespaces.GetFESpaceAtLevel(level), level);
      }
   }

   virtual ~ElasticityMultigrid()
   {
      delete amg;
   }

private:
   void ConstructBilinearForm(ParFiniteElementSpace& fespace,
                              bool partial_assembly)
   {
      ParBilinearForm* form = new ParBilinearForm(&fespace);
      if (partial_assembly)
      {
         form->SetAssemblyLevel(AssemblyLevel::PARTIAL);
      }
      if(use_optimized_kernel)
      {
         form->AddDomainIntegrator(new MyElasticityIntegrator(lambda_func, mu_func));
      }
      else
      {
         form->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));
      }
      form->Assemble();
      bfs.Append(form);
   }

   void ConstructCoarseOperatorAndSolver(ParFiniteElementSpace& coarse_fespace)
   {
      ConstructBilinearForm(coarse_fespace, false);

      HypreParMatrix* hypreCoarseMat = new HypreParMatrix();
      bfs[0]->FormSystemMatrix(*essentialTrueDofs[0], *hypreCoarseMat);

      amg = new HypreBoomerAMG(*hypreCoarseMat);
      amg->SetPrintLevel(-1);
      int vdim = coarse_fespace.GetVDim();
      // amg->SetSystemsOptions(vdim);

      // amg->SetElasticityOptions(&coarse_fespace);
      amg->SetSystemsOptions(vdim, true);

      // TODO: Strong threshold
      // amg->SetStrongThreshold(0.9);

      CGSolver* pcg = new CGSolver(MPI_COMM_WORLD);
      pcg->SetPrintLevel(-1);
      pcg->SetMaxIter(10);
      pcg->SetRelTol(sqrt(1e-4));
      pcg->SetAbsTol(0.0);
      pcg->SetOperator(*hypreCoarseMat);
      pcg->SetPreconditioner(*amg);

      AddLevel(hypreCoarseMat, pcg, true, true);
   }

   void ConstructOperatorAndSmoother(ParFiniteElementSpace& fespace, int level)
   {
      const Array<int> &ess_tdof_list = *essentialTrueDofs[level];
      ConstructBilinearForm(fespace, pa);
      if(pa)
      {
         OperatorPtr opr;
         opr.SetType(Operator::ANY_TYPE);
         bfs.Last()->FormSystemMatrix(ess_tdof_list, opr);
         ownOperator = opr.OwnsOperator();
         opr.SetOperatorOwner(false);

         Vector diag(fespace.GetTrueVSize());
         bfs.Last()->AssembleDiagonal(diag);

         Solver* smoother = new OperatorChebyshevSmoother(
            *opr, diag, ess_tdof_list, 2, fespace.GetParMesh()->GetComm());

         AddLevel(opr.Ptr(), smoother, ownOperator, true);
      }
      else
      {
         HypreParMatrix* hypreMat = new HypreParMatrix();
         bfs.Last()->FormSystemMatrix(ess_tdof_list, *hypreMat);

         Vector diag(fespace.GetTrueVSize());
         bfs.Last()->AssembleDiagonal(diag);

         Solver* smoother = new OperatorChebyshevSmoother(
             *hypreMat, diag, ess_tdof_list, 2, fespace.GetParMesh()->GetComm());

         AddLevel(hypreMat, smoother, true, true);
      }
   }
};

int main(int argc, char *argv[])
{
   // Initialize MPI and HYPRE
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   MPI_Barrier(MPI_COMM_WORLD);

   StopWatch total_timer;
   total_timer.Start();

   StopWatch setup_timer;
   setup_timer.Start();

   // Parse options
   const char *mesh_file = "../mfem/data/beam-hex.mesh";
   int order = 1;
   int par_ref_levels = 2;
   int geometric_refinements = 0;
   int order_refinements = 0;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = false;
   bool op = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&par_ref_levels, "-r", "--refine-parallel",
                  "Number of uniform parallel refinements.");
   args.AddOption(&geometric_refinements, "-gr", "--geometric-refinements",
                  "Number of geometric refinements done prior to order refinements.");
   args.AddOption(&order_refinements, "-or", "--order-refinements",
                  "Finest polynomial order in the multigrid hierarchy.");
   args.AddOption(&op, "-op", "--optimize", "-no-op",
                  "--no-optimize", "Enable the optimized elasticity integrator.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // Configure device
   Device device(device_config);
   if (myid == 0)
   {
      device.Print();
      if(pa)
      {
         cout << "Partial assembly" << endl;
      }
      else
      {
         cout << "Full assembly" << endl;
      }
   }

   // Read mesh
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   MFEM_VERIFY(mesh->SpaceDimension() == dim, "invalid mesh");

   if (mesh->attributes.Max() < 2 || mesh->bdr_attributes.Max() < 2)
   {
      cerr << "\nInput mesh should have at least two materials and "
           << "two boundary attributes! (See schematic in ex2.cpp)\n"
           << endl;
      return 3;
   }

   // Check affine geometry
   if (pa && op && (geometric_refinements > 0 || order_refinements > 1))
   {
      setup_timer.Stop();
      total_timer.Stop();
      if (myid == 0)
      {
         checkAffineMesh(*mesh);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      total_timer.Start();
      setup_timer.Start();
   }

   // Refine mesh
   {
      int ref_levels =
         (int)floor(log(1000./mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }
   mesh->EnsureNCMesh();

   // Partition mesh
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   {
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
   }

   // Create finite element space
   FiniteElementCollection *fec = new H1_FECollection(order, dim);
   ParFiniteElementSpace *coarse_fespace = new ParFiniteElementSpace(pmesh,
                fec, dim, Ordering::byNODES);

   Array<FiniteElementCollection*> collections;
   collections.Append(fec);

   ParFiniteElementSpaceHierarchy* fespaces = new ParFiniteElementSpaceHierarchy(
      pmesh, coarse_fespace, true, true);

   // Geometric refinement
   for (int level = 0; level < geometric_refinements; ++level)
   {
      fespaces->AddUniformlyRefinedLevel(dim, Ordering::byNODES);
   }

   // Polynomial hierarchy
   const int finest_order = order_refinements;

   std::vector<int> p_hierarchy;
   int p = finest_order;
   while (p > 0)
   {
      if (p == 1) { break; }
      p_hierarchy.push_back(p);
      int next_p = std::max(1, p / 2);
      if (next_p >= p) { next_p = p - 1; }
      p = next_p;
   }

   std::reverse(p_hierarchy.begin(), p_hierarchy.end());

   // Order refinement
   for (int refined_order : p_hierarchy)
   {
      collections.Append(new H1_FECollection(refined_order, dim));
      fespaces->AddOrderRefinedLevel(collections.Last(), dim, Ordering::byNODES);
   }

   // Print hierarchy
   HYPRE_BigInt size = fespaces->GetFinestFESpace().GlobalTrueVSize();

   if (myid == 0)
   {
      cout << "=========================" << endl;
      cout << "Multigrid hierarchy:" << endl;
      cout << "Levels: " << fespaces->GetNumLevels() << endl;
      cout << "Total DOFs: " << size << endl;
      for (int level = 0; level < fespaces->GetNumLevels(); ++level)
      {
         ParFiniteElementSpace &space = fespaces->GetFESpaceAtLevel(level);
         HYPRE_BigInt level_dofs = space.GlobalTrueVSize();
         int level_order = space.FEColl()->GetOrder();
         int num_elements = space.GetNE();
         cout << "Level " << level << ": order = " << level_order
              << ", global DOFs = " << level_dofs
              << ", local elements = " << num_elements << endl;
      }
      cout << "=========================" << endl;
   }

   // Boundary load
   VectorArrayCoefficient f(dim);
   for (int i = 0; i < dim-1; i++)
   {
      f.Set(i, new ConstantCoefficient(0.0));
   }
   {
      Vector pull_force(pmesh->bdr_attributes.Max());
      pull_force = 0.0;
      pull_force(1) = -1.0e-2;
      f.Set(dim-1, new PWConstCoefficient(pull_force));
   }

   ParLinearForm *b = new ParLinearForm(&fespaces->GetFinestFESpace());
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));

   // Material coefficients
   Vector lambda(pmesh->attributes.Max());
   lambda = 1.0;
   lambda(0) = lambda(1)*50;
   PWConstCoefficient lambda_func(lambda);
   Vector mu(pmesh->attributes.Max());
   mu = 1.0;
   mu(0) = mu(1)*50;
   PWConstCoefficient mu_func(mu);

   if(op)
   {
      BilinearFormIntegrator *integ = new MyElasticityIntegrator(lambda_func,mu_func);
   }
   else
   {
      BilinearFormIntegrator *integ = new ElasticityIntegrator(lambda_func,mu_func);
   }

   // Initialize solution
   Vector zero_vec(dim);
   zero_vec = 0.0;
   VectorConstantCoefficient zero_vec_coeff(zero_vec);
   ParGridFunction x(&fespaces->GetFinestFESpace());
   x = 0.0;

   // Essential boundaries
   Array<int> ess_bdr(pmesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;

   setup_timer.Stop();

   StopWatch assemble_timer;
   StopWatch setup_muitgrid_timer;
   StopWatch formfls_timer;
   StopWatch solver_timer;

   assemble_timer.Start();

   // Assemble right-hand side
   b->Assemble();

   assemble_timer.Stop();

   // Apply boundary values
   x.ProjectBdrCoefficient(zero_vec_coeff, ess_bdr);

   setup_muitgrid_timer.Start();

   // Construct multigrid
   ElasticityMultigrid* M = new ElasticityMultigrid(*fespaces, ess_bdr, lambda_func, mu_func, pa, op);
   M->SetCycleType(Multigrid::CycleType::VCYCLE, 1, 1);

   setup_muitgrid_timer.Stop();

   OperatorPtr A;
   Vector B, X;

   formfls_timer.Start();

   // Form linear system
   M->FormFineLinearSystem(x, *b, A, X, B);

   formfls_timer.Stop();

   solver_timer.Start();

   // Solve
   CGSolver pcg(MPI_COMM_WORLD);
   pcg.SetRelTol(1e-6);
   pcg.SetMaxIter(2000);
   pcg.SetPrintLevel(1);
   pcg.SetOperator(*A);
   pcg.SetPreconditioner(*M);

   pcg.Mult(B, X);

   solver_timer.Stop();

   if (myid == 0)
   {
      cout << "Setup time: " << setup_timer.RealTime() << " seconds" << endl;
      cout << "Right-hand side assembly time: " << assemble_timer.RealTime()
           << " seconds" << endl;
      cout << "Multigrid setup time: " << setup_muitgrid_timer.RealTime()
           << " seconds" << endl;
      cout << "Linear system assembly time: " << formfls_timer.RealTime()
           << " seconds" << endl;
      cout << "Solve time: " << solver_timer.RealTime() << " seconds" << endl;
      cout << "PCG iterations: " << pcg.GetNumIterations() << endl;
      cout << "PCG relative residual: " << pcg.GetFinalRelNorm() << endl;
      cout << "Combined time: " << setup_timer.RealTime() + assemble_timer.RealTime() +
              setup_muitgrid_timer.RealTime() + formfls_timer.RealTime() +
              solver_timer.RealTime() << " seconds" << endl;
   }

   // Recover solution
   M->RecoverFineFEMSolution(X, *b, x);

   // Save mesh and solution
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      fespaces->GetFinestFESpace().GetParMesh()->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // Visualize
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *fespaces->GetFinestFESpace().GetParMesh()
               << x << flush;
   }

   total_timer.Stop();
   if (myid == 0)
      cout << "Total time: " << total_timer.RealTime() << " seconds" << endl;

   // Clean up
   delete M;
   delete b;
   delete fespaces;
   for (int level = 0; level < collections.Size(); ++level)
   {
      delete collections[level];
   }

   return 0;
}
