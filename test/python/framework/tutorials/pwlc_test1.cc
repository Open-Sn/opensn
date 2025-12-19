// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "test/python/src/bindings.h"

using namespace opensn;

namespace unit_tests
{

/**
 * This is a simple test of the Finite Volume spatial discretization applied
 * to Laplace's problem.
 */

void
SimTest03_PWLC(std::shared_ptr<MeshContinuum> grid)
{
  opensn::log.Log() << "Coding Tutorial 3";

  opensn::log.Log() << "Global num cells: " << grid->GetGlobalNumberOfCells();

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearContinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const auto num_global_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_global_dofs;

  // Initializes Mats and Vecs
  const auto n = static_cast<PetscInt>(num_local_dofs);
  const auto N = static_cast<PetscInt>(num_global_dofs);
  Mat A;
  Vec x, b;

  A = CreateSquareMatrix(n, N);
  x = CreateVector(n, N);
  b = CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  InitMatrixSparsity(A, nodal_nnz_in_diag, nodal_nnz_off_diag);

  // Assemble the system
  opensn::log.Log() << "Assembling system: ";
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    const size_t num_nodes = cell_mapping.GetNumNodes();
    DenseMatrix<double> Acell(num_nodes, num_nodes, 0.0);
    Vector<double> cell_rhs(num_nodes, 0.0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
      {
        double entry_aij = 0.0;
        for (size_t qp : fe_vol_data.GetQuadraturePointIndices())
        {
          entry_aij +=
            fe_vol_data.ShapeGrad(i, qp).Dot(fe_vol_data.ShapeGrad(j, qp)) * fe_vol_data.JxW(qp);
        } // for qp
        Acell(i, j) = entry_aij;
      } // for j
      for (size_t qp : fe_vol_data.GetQuadraturePointIndices())
        cell_rhs(i) += 1.0 * fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
    } // for i

    // Flag nodes for being on dirichlet boundary
    std::vector<bool> node_boundary_flag(num_nodes, false);
    const size_t num_faces = cell.faces.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      if (face.has_neighbor)
        continue;

      const size_t num_face_nodes = face.vertex_ids.size();
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const uint i = cell_mapping.MapFaceNode(f, fi);
        node_boundary_flag[i] = true;
      } // for fi
    } // for face f

    // Develop node mapping
    std::vector<uint64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    // Assembly into system
    for (size_t i = 0; i < num_nodes; ++i)
    {
      if (node_boundary_flag[i]) // if dirichlet node
      {
        MatSetValue(A, imap[i], imap[i], 1.0, ADD_VALUES);
        VecSetValue(b, imap[i], 0.0, ADD_VALUES);
      }
      else
      {
        for (size_t j = 0; j < num_nodes; ++j)
        {
          if (not node_boundary_flag[j])
            MatSetValue(A, imap[i], imap[j], Acell(i, j), ADD_VALUES);
        } // for j
        VecSetValue(b, imap[i], cell_rhs(i), ADD_VALUES);
      }
    } // for i
  } // for cell

  opensn::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  opensn::log.Log() << "Done global assembly";

  // Create Krylov Solver
  opensn::log.Log() << "Solving: ";
  auto petsc_solver =
    CreateCommonKrylovSolverSetup(A, "PWLCDiffSolver", KSPCG, PCGAMG, 1.0e-15, 0.0, 30);

  // Solve
  KSPSolve(petsc_solver.ksp, b, x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(petsc_solver.ksp, &reason);
  if (reason == KSP_CONVERGED_RTOL)
    opensn::log.Log() << "Converged";

  opensn::log.Log() << "Done solving";

  // Extract PETSc vector
  std::vector<double> field;
  sdm.LocalizePETScVector(x, field, OneDofPerNode);

  // Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  opensn::log.Log() << "Done cleanup";

  // Create Field Function
  auto ff = std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, Unknown(UnknownType::SCALAR));

  ff->UpdateFieldVector(field);

  FieldFunctionGridBased::ExportMultipleToPVTU("CodeTut3_PWLC", {ff});
}

} // namespace unit_tests
