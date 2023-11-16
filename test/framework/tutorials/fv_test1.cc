#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/petsc_utils/petsc_utils.h"

#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/physics/field_function/field_function_grid_based.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/console/console.h"

using namespace opensn;

namespace chi_unit_sim_tests
{

/**This is a simple test of the Finite Volume spatial discretization applied
 * to Laplace's problem. */
ParameterBlock chiSimTest01_FV(const InputParameters& params);

RegisterWrapperFunction(chi_unit_sim_tests, chiSimTest01_FV, nullptr, chiSimTest01_FV);

ParameterBlock
chiSimTest01_FV(const InputParameters&)
{
  opensn::Chi::log.Log() << "Coding Tutorial 1";

  // Get grid
  auto grid_ptr = GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  opensn::Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make SDM
  typedef std::shared_ptr<SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = FiniteVolume::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::Chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  opensn::Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

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
  opensn::Chi::log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const int64_t imap = sdm.MapDOF(cell, 0);

    const auto& xp = cell.centroid_;
    const double V = cell_mapping.CellVolume();

    size_t f = 0;
    for (const auto& face : cell.faces_)
    {
      const auto Af = face.normal_ * cell_mapping.FaceArea(f);

      if (face.has_neighbor_)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id_];
        const int64_t jnmap = sdm.MapDOF(adj_cell, 0);

        const auto& xn = adj_cell.centroid_;

        const auto xpn = xn - xp;

        const auto cf = Af.Dot(xpn) / xpn.NormSquare();

        MatSetValue(A, imap, imap, cf, ADD_VALUES);
        MatSetValue(A, imap, jnmap, -cf, ADD_VALUES);
      }
      else
      {
        const auto& xn = xp + 2.0 * (face.centroid_ - xp);
        const auto xpn = xn - xp;

        const auto cf = Af.Dot(xpn) / xpn.NormSquare();

        MatSetValue(A, imap, imap, cf, ADD_VALUES);
      }
      ++f;
    } // for face

    VecSetValue(b, imap, 1.0 * V, ADD_VALUES);
  } // for cell i

  opensn::Chi::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  opensn::Chi::log.Log() << "Done global assembly";

  // Create Krylov Solver
  opensn::Chi::log.Log() << "Solving: ";
  auto petsc_solver = CreateCommonKrylovSolverSetup(A, "FVDiffSolver", KSPCG, PCGAMG, 1.0e-6, 1000);

  // Solve
  KSPSolve(petsc_solver.ksp, b, x);

  opensn::Chi::log.Log() << "Done solving";

  // Extract PETSc vector
  std::vector<double> field(num_local_dofs, 0.0);
  sdm.LocalizePETScVector(x, field, OneDofPerNode);

  // Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  opensn::Chi::log.Log() << "Done cleanup";

  // Create Field Function
  auto ff = std::make_shared<FieldFunctionGridBased>("Phi", sdm_ptr, Unknown(UnknownType::SCALAR));

  ff->UpdateFieldVector(field);

  FieldFunctionGridBased::ExportMultipleToVTK("CodeTut1_FV", {ff});

  return ParameterBlock();
}

} // namespace chi_unit_sim_tests
