// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/vector_ghost_communicator/vector_ghost_communicator.h"
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
SimTest02_FV(std::shared_ptr<MeshContinuum> grid)
{
  opensn::log.Log() << "Coding Tutorial 2";

  opensn::log.Log() << "Global num cells: " << grid->GetGlobalNumberOfCells();

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = FiniteVolume::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_global_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_global_dofs;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_global_dofs);
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
    const int64_t imap = sdm.MapDOF(cell, 0);

    const auto& xp = cell.centroid;
    const auto V = cell.volume;

    for (const auto& face : cell.faces)
    {
      const auto Af = face.normal * face.area;

      if (face.has_neighbor)
      {
        const auto& adj_cell = grid->cells[face.neighbor_id];
        const int64_t jnmap = sdm.MapDOF(adj_cell, 0);

        const auto& xn = adj_cell.centroid;

        const auto xpn = xn - xp;

        const auto cf = Af.Dot(xpn) / xpn.NormSquare();

        MatSetValue(A, imap, imap, cf, ADD_VALUES);
        MatSetValue(A, imap, jnmap, -cf, ADD_VALUES);
      }
      else
      {
        const auto& xn = xp + 2.0 * (face.centroid - xp);
        const auto xpn = xn - xp;

        const auto cf = Af.Dot(xpn) / xpn.NormSquare();

        MatSetValue(A, imap, imap, cf, ADD_VALUES);
      }
    } // for face

    VecSetValue(b, imap, 1.0 * V, ADD_VALUES);
  } // for cell i

  opensn::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  opensn::log.Log() << "Done global assembly";

  // Create Krylov Solver
  opensn::log.Log() << "Solving: ";
  auto petsc_solver =
    CreateCommonKrylovSolverSetup(A, "FVDiffSolver", KSPCG, PCGAMG, 0.0, 1.0e-9, 18);

  // Solve
  KSPSolve(petsc_solver.ksp, b, x);
  KSPConvergedReason reason;
  KSPGetConvergedReason(petsc_solver.ksp, &reason);
  if (reason == KSP_CONVERGED_ATOL)
    opensn::log.Log() << "Converged";

  // Extract PETSc vector
  std::vector<double> field(num_local_dofs, 0.0);
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

  FieldFunctionGridBased::ExportMultipleToVTK("CodeTut2_FV", {ff});

  // Make ghosted vectors
  std::vector<int64_t> ghost_ids = sdm.GetGhostDOFIndices(OneDofPerNode);

  VectorGhostCommunicator vgc(num_local_dofs, num_global_dofs, ghost_ids, opensn::mpi_comm);
  std::vector<double> field_wg = vgc.MakeGhostedVector(field);

  vgc.CommunicateGhostEntries(field_wg);

  // Setup gradient unknown
  // structure
  UnknownManager grad_uk_man({Unknown(UnknownType::VECTOR_3)});

  const size_t num_grad_dofs = sdm.GetNumLocalDOFs(grad_uk_man);

  std::vector<double> grad_phi(num_grad_dofs, 0.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);

    const int64_t imap = sdm.MapDOFLocal(cell, 0);
    const double phi_P = field_wg[imap];

    const auto& xp = cell.centroid;

    auto grad_phi_P = Vector3(0, 0, 0);

    for (const auto& face : cell.faces)
    {
      const auto& xf = face.centroid;
      const auto Af = face.normal * face.area;

      double phi_N = 0.0;
      auto xn = xp + 2 * (xf - xp);

      if (face.has_neighbor)
      {
        const auto& adj_cell = grid->cells[face.neighbor_id];
        const int64_t nmap = sdm.MapDOFLocal(adj_cell, 0);
        phi_N = field_wg[nmap];

        xn = adj_cell.centroid;
      }

      grad_phi_P += Af * ((xn - xf).Norm() * phi_P + (xf - xp).Norm() * phi_N) / (xn - xp).Norm();
    } // for face
    grad_phi_P /= cell.volume;

    const int64_t xmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 0);
    const int64_t ymap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 1);
    const int64_t zmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 2);

    grad_phi[xmap] = grad_phi_P.x;
    grad_phi[ymap] = grad_phi_P.y;
    grad_phi[zmap] = grad_phi_P.z;
  } // for cell

  // Create Field Function
  auto ff_grad =
    std::make_shared<FieldFunctionGridBased>("GradPhi", sdm_ptr, Unknown(UnknownType::VECTOR_3));

  ff_grad->UpdateFieldVector(grad_phi);

  FieldFunctionGridBased::ExportMultipleToVTK("CodeTut2_FV_grad", {ff_grad});
}

} // namespace unit_tests
