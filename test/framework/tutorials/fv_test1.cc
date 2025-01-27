#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/lib/console.h"

using namespace opensn;

namespace unit_sim_tests
{

/**
 * This is a simple test of the Finite Volume spatial discretization applied
 * to Laplace's problem.
 */

void
SimTest01_FV()
{
  opensn::log.Log() << "Coding Tutorial 1";

  // Get grid
  auto grid_ptr = GetCurrentMesh();
  const auto& grid = *grid_ptr;

  opensn::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = FiniteVolume::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);
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
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const int64_t imap = sdm.MapDOF(cell, 0);

    const auto& xp = cell.centroid;
    const double V = cell_mapping.GetCellVolume();

    size_t f = 0;
    for (const auto& face : cell.faces)
    {
      const auto Af = face.normal * cell_mapping.GetFaceArea(f);

      if (face.has_neighbor)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id];
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
      ++f;
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
    CreateCommonKrylovSolverSetup(A, "FVDiffSolver", KSPCG, PCGAMG, 0.0, 1.0e-15, 25);

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

  FieldFunctionGridBased::ExportMultipleToVTK("CodeTut1_FV", {ff});
}

BIND_FUNCTION(unit_sim_tests, SimTest01_FV);

} // namespace unit_sim_tests
