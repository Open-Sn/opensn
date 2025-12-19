// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/diffusion/diffusion.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

DiffusionSolver::DiffusionSolver(std::string name,
                                 const opensn::SpatialDiscretization& sdm,
                                 const UnknownManager& uk_man,
                                 std::map<uint64_t, BoundaryCondition> bcs,
                                 MatID2XSMap map_mat_id_2_xs,
                                 const std::vector<UnitCellMatrices>& unit_cell_matrices,
                                 const bool suppress_bcs,
                                 const bool requires_ghosts,
                                 const bool verbose)
  : name_(std::move(name)),
    grid_(sdm.GetGrid()),
    sdm_(sdm),
    uk_man_(uk_man),
    bcs_(std::move(bcs)),
    mat_id_2_xs_map_(std::move(map_mat_id_2_xs)),
    unit_cell_matrices_(unit_cell_matrices),
    num_local_dofs_(static_cast<PetscInt>(sdm_.GetNumLocalDOFs(uk_man_))),
    num_global_dofs_(static_cast<PetscInt>(sdm_.GetNumGlobalDOFs(uk_man_))),
    A_(nullptr),
    rhs_(nullptr),
    ksp_(nullptr),
    requires_ghosts_(requires_ghosts),
    suppress_bcs_(suppress_bcs)
{
  options.verbose = verbose;
}

DiffusionSolver::~DiffusionSolver()
{
  OpenSnPETScCall(MatDestroy(&A_));
  OpenSnPETScCall(VecDestroy(&rhs_));
  OpenSnPETScCall(KSPDestroy(&ksp_));
}

std::string
DiffusionSolver::GetName() const
{
  return name_;
}

const Vec&
DiffusionSolver::GetRHS() const
{
  return rhs_;
}

const UnknownManager&
DiffusionSolver::GetUnknownStructure() const
{
  return uk_man_;
}

const SpatialDiscretization&
DiffusionSolver::GetSpatialDiscretization() const
{
  return sdm_;
}

std::pair<std::uint64_t, std::uint64_t>
DiffusionSolver::GetNumPhiIterativeUnknowns()
{
  return {sdm_.GetNumLocalDOFs(uk_man_), sdm_.GetNumGlobalDOFs(uk_man_)};
}

void
DiffusionSolver::AddToRHS(const std::vector<double>& values)
{
  const auto num_local_dofs = sdm_.GetNumLocalDOFs(uk_man_);
  if (num_local_dofs != values.size())
    throw std::invalid_argument("Vector size mismatch.");

  PetscScalar* rhs_ptr = nullptr;
  OpenSnPETScCall(VecGetArray(rhs_, &rhs_ptr));
  for (std::uint64_t i = 0; i < num_local_dofs; ++i)
    rhs_ptr[i] += values[i];
  OpenSnPETScCall(VecRestoreArray(rhs_, &rhs_ptr));
}

void
DiffusionSolver::AddToMatrix(const std::vector<PetscInt>& rows,
                             const std::vector<PetscInt>& cols,
                             const std::vector<double>& vals)
{
  if (rows.size() != cols.size() or rows.size() != vals.size())
    throw std::invalid_argument("The number of row entries, column entries, and value "
                                "entries do not agree.");
  for (int i = 0; i < vals.size(); ++i)
    OpenSnPETScCall(MatSetValue(A_, rows[i], cols[i], vals[i], ADD_VALUES));
  OpenSnPETScCall(MatAssemblyBegin(A_, MAT_FLUSH_ASSEMBLY));
  OpenSnPETScCall(MatAssemblyEnd(A_, MAT_FLUSH_ASSEMBLY));
}

void
DiffusionSolver::Initialize()
{
  if (options.verbose)
    log.Log() << name_ << ": Initializing PETSc items";

  if (options.verbose)
    log.Log() << name_ << ": Global number of DOFs=" << num_global_dofs_;

  opensn::mpi_comm.barrier();
  log.Log() << "Sparsity pattern";
  opensn::mpi_comm.barrier();
  // Create Matrix
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm_.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, uk_man_);
  opensn::mpi_comm.barrier();
  log.Log() << "Done Sparsity pattern";
  opensn::mpi_comm.barrier();
  A_ = CreateSquareMatrix(num_local_dofs_, num_global_dofs_);
  InitMatrixSparsity(A_, nodal_nnz_in_diag, nodal_nnz_off_diag);
  opensn::mpi_comm.barrier();
  log.Log() << "Done matrix creation";
  opensn::mpi_comm.barrier();

  // Create RHS
  if (not requires_ghosts_)
    rhs_ = CreateVector(num_local_dofs_, num_global_dofs_);
  else
  {
    auto ghost_ids = sdm_.GetGhostDOFIndices(uk_man_);
    std::vector<int64_t> ghids(ghost_ids.begin(), ghost_ids.end());
    rhs_ = CreateVectorWithGhosts(num_local_dofs_,
                                  num_global_dofs_,
                                  static_cast<PetscInt>(sdm_.GetNumGhostDOFs(uk_man_)),
                                  ghids);
  }

  opensn::mpi_comm.barrier();
  log.Log() << "Done vector creation";
  opensn::mpi_comm.barrier();

  // Create KSP
  OpenSnPETScCall(KSPCreate(opensn::mpi_comm, &ksp_));
  OpenSnPETScCall(KSPSetOptionsPrefix(ksp_, name_.c_str()));
  OpenSnPETScCall(KSPSetType(ksp_, KSPCG));

  OpenSnPETScCall(
    KSPSetTolerances(ksp_, 1.0e-50, options.residual_tolerance, 1.0e50, options.max_iters));

  // Set Pre-conditioner
  PC pc = nullptr;
  OpenSnPETScCall(KSPGetPC(ksp_, &pc));
  //  PCSetType(pc, PCGAMG);
  OpenSnPETScCall(PCSetType(pc, PCHYPRE));

  OpenSnPETScCall(PCHYPRESetType(pc, "boomeramg"));
  std::vector<std::string> pc_options = {"pc_hypre_boomeramg_agg_nl 1",
                                         "pc_hypre_boomeramg_P_max 4",
                                         "pc_hypre_boomeramg_grid_sweeps_coarse 1",
                                         "pc_hypre_boomeramg_max_levels 25",
                                         "pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi",
                                         "pc_hypre_boomeramg_coarsen_type HMIS",
                                         "pc_hypre_boomeramg_interp_type ext+i"};

  if (grid_->GetDimension() == 2)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
  else if (grid_->GetDimension() == 3)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.8");

  for (const auto& option : pc_options)
    OpenSnPETScCall(PetscOptionsInsertString(nullptr, ("-" + name_ + option).c_str()));

  OpenSnPETScCall(PetscOptionsInsertString(nullptr, options.additional_options_string.c_str()));

  OpenSnPETScCall(PCSetFromOptions(pc));
  OpenSnPETScCall(KSPSetFromOptions(ksp_));
}

void
DiffusionSolver::Solve(std::vector<double>& solution, bool use_initial_guess)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::Solve";
  Vec x = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &x));
  OpenSnPETScCall(VecSet(x, 0.0));

  if (not use_initial_guess)
    OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE));
  else
    OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE));

  OpenSnPETScCall(KSPSetTolerances(
    ksp_, options.residual_tolerance, options.residual_tolerance, 1.0e50, options.max_iters));

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    OpenSnPETScCall(MatIsSymmetric(A_, 1.0e-6, &symmetry));
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    OpenSnPETScCall(KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr));

    double rhs_norm = 0.0;
    OpenSnPETScCall(VecNorm(rhs_, NORM_2, &rhs_norm));
    log.Log() << "RHS-norm " << rhs_norm;
  }

  if (use_initial_guess)
  {
    double* x_raw = nullptr;
    OpenSnPETScCall(VecGetArray(x, &x_raw));
    size_t k = 0;
    for (const auto& value : solution)
      x_raw[k++] = value;
    OpenSnPETScCall(VecRestoreArray(x, &x_raw));
  }

  // Solve
  OpenSnPETScCall(KSPSolve(ksp_, rhs_, x));

  // Print convergence info
  if (options.verbose)
  {
    double sol_norm = 0.0;
    OpenSnPETScCall(VecNorm(x, NORM_2, &sol_norm));
    log.Log() << "Solution-norm " << sol_norm;

    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    OpenSnPETScCall(KSPGetConvergedReason(ksp_, &reason));

    log.Log() << "Convergence Reason: " << GetPETScConvergedReasonstring(reason);
  }

  // Transfer petsc solution to vector
  if (requires_ghosts_)
  {
    CommunicateGhostEntries(x);
    sdm_.LocalizePETScVectorWithGhosts(x, solution, uk_man_);
  }
  else
    sdm_.LocalizePETScVector(x, solution, uk_man_);

  // Cleanup x
  OpenSnPETScCall(VecDestroy(&x));
}

void
DiffusionSolver::Solve(Vec petsc_solution, bool use_initial_guess)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::Solve";
  Vec x = nullptr;
  OpenSnPETScCall(VecDuplicate(rhs_, &x));
  OpenSnPETScCall(VecSet(x, 0.0));

  if (not use_initial_guess)
    OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE));
  else
    OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE));

  OpenSnPETScCall(KSPSetTolerances(
    ksp_, options.residual_tolerance, options.residual_tolerance, 1.0e50, options.max_iters));

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    OpenSnPETScCall(MatIsSymmetric(A_, 1.0e-6, &symmetry));
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    OpenSnPETScCall(KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr));

    double rhs_norm = 0.0;
    OpenSnPETScCall(VecNorm(rhs_, NORM_2, &rhs_norm));
    log.Log() << "RHS-norm " << rhs_norm;
  }

  if (use_initial_guess)
  {
    OpenSnPETScCall(VecCopy(petsc_solution, x));
  }

  // Solve
  OpenSnPETScCall(KSPSolve(ksp_, rhs_, x));

  // Print convergence info
  if (options.verbose)
  {
    double sol_norm = 0.0;
    OpenSnPETScCall(VecNorm(x, NORM_2, &sol_norm));
    log.Log() << "Solution-norm " << sol_norm;

    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    OpenSnPETScCall(KSPGetConvergedReason(ksp_, &reason));

    log.Log() << "Convergence Reason: " << GetPETScConvergedReasonstring(reason);
  }

  // Transfer petsc solution to vector
  OpenSnPETScCall(VecCopy(x, petsc_solution));

  // Cleanup x
  OpenSnPETScCall(VecDestroy(&x));
}

} // namespace opensn
