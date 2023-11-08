#include "framework/mesh/MeshHandler/chi_meshhandler.h"

#include "framework/math/SpatialDiscretization/FiniteVolume/FiniteVolume.h"
#include "framework/math/PETScUtils/petsc_utils.h"

#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/physics/FieldFunction/fieldfunction_gridbased.h"

#include "framework/math/VectorGhostCommunicator/vector_ghost_communicator.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/console/chi_console.h"

namespace chi_unit_sim_tests
{

chi::ParameterBlock chiSimTest02_FV(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_sim_tests,
                        /*name_in_lua=*/chiSimTest02_FV,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chiSimTest02_FV);

/**This is a simple test of the Finite Volume spatial discretization applied
 * to Laplace's problem. */
chi::ParameterBlock
chiSimTest02_FV(const chi::InputParameters&)
{
  Chi::log.Log() << "Coding Tutorial 2";

  //============================================= Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= Make SDM
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = chi_math::spatial_discretization::FiniteVolume::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  Chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);
  Mat A;
  Vec x, b;

  A = chi_math::PETScUtils::CreateSquareMatrix(n, N);
  x = chi_math::PETScUtils::CreateVector(n, N);
  b = chi_math::PETScUtils::CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  chi_math::PETScUtils::InitMatrixSparsity(A, nodal_nnz_in_diag, nodal_nnz_off_diag);

  //============================================= Assemble the system
  Chi::log.Log() << "Assembling system: ";
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

  Chi::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  Chi::log.Log() << "Done global assembly";

  //============================================= Create Krylov Solver
  Chi::log.Log() << "Solving: ";
  auto petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(A,              // Matrix
                                                        "FVDiffSolver", // Solver name
                                                        KSPCG,          // Solver type
                                                        PCGAMG,         // Preconditioner type
                                                        1.0e-6, // Relative residual tolerance
                                                        1000);  // Max iterations

  //============================================= Solve
  KSPSolve(petsc_solver.ksp, b, x);

  Chi::log.Log() << "Done solving";

  //============================================= Extract PETSc vector
  std::vector<double> field(num_local_dofs, 0.0);
  sdm.LocalizePETScVector(x, field, OneDofPerNode);

  //============================================= Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  Chi::log.Log() << "Done cleanup";

  //============================================= Create Field Function
  auto ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
    "Phi", sdm_ptr, chi_math::Unknown(chi_math::UnknownType::SCALAR));

  ff->UpdateFieldVector(field);

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK("CodeTut2_FV", {ff});

  //============================================= Make ghosted vectors
  std::vector<int64_t> ghost_ids = sdm.GetGhostDOFIndices(OneDofPerNode);

  chi_math::VectorGhostCommunicator vgc(num_local_dofs, num_globl_dofs, ghost_ids, Chi::mpi.comm);
  std::vector<double> field_wg = vgc.MakeGhostedVector(field);

  vgc.CommunicateGhostEntries(field_wg);

  //============================================= Setup gradient unknown
  // structure
  chi_math::UnknownManager grad_uk_man({chi_math::Unknown(chi_math::UnknownType::VECTOR_3)});

  const size_t num_grad_dofs = sdm.GetNumLocalDOFs(grad_uk_man);

  std::vector<double> grad_phi(num_grad_dofs, 0.0);

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);

    const int64_t imap = sdm.MapDOFLocal(cell, 0);
    const double phi_P = field_wg[imap];

    const auto& xp = cell.centroid_;

    auto grad_phi_P = chi_mesh::Vector3(0, 0, 0);

    size_t f = 0;
    for (const auto& face : cell.faces_)
    {
      const auto& xf = face.centroid_;
      const auto Af = cell_mapping.FaceArea(f) * face.normal_;

      double phi_N = 0.0;
      auto xn = xp + 2 * (xf - xp);

      if (face.has_neighbor_)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id_];
        const int64_t nmap = sdm.MapDOFLocal(adj_cell, 0);
        phi_N = field_wg[nmap];

        xn = adj_cell.centroid_;
      }

      grad_phi_P += Af * ((xn - xf).Norm() * phi_P + (xf - xp).Norm() * phi_N) / (xn - xp).Norm();
      ++f;
    } // for face
    grad_phi_P /= cell_mapping.CellVolume();

    const int64_t xmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 0);
    const int64_t ymap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 1);
    const int64_t zmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 2);

    grad_phi[xmap] = grad_phi_P.x;
    grad_phi[ymap] = grad_phi_P.y;
    grad_phi[zmap] = grad_phi_P.z;
  } // for cell

  //============================================= Create Field Function
  auto ff_grad = std::make_shared<chi_physics::FieldFunctionGridBased>(
    "GradPhi", sdm_ptr, chi_math::Unknown(chi_math::UnknownType::VECTOR_3));

  ff_grad->UpdateFieldVector(grad_phi);

  chi_physics::FieldFunctionGridBased::ExportMultipleToVTK("CodeTut2_FV_grad", {ff_grad});

  return chi::ParameterBlock();
}

} // namespace chi_unit_sim_tests
