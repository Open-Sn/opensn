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
    num_local_dofs_(static_cast<int64_t>(sdm_.GetNumLocalDOFs(uk_man_))),
    num_global_dofs_(static_cast<int64_t>(sdm_.GetNumGlobalDOFs(uk_man_))),
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
  MatDestroy(&A_);
  VecDestroy(&rhs_);
  KSPDestroy(&ksp_);
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

std::pair<size_t, size_t>
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

  PetscScalar* rhs_ptr;
  VecGetArray(rhs_, &rhs_ptr);
  for (size_t i = 0; i < num_local_dofs; ++i)
    rhs_ptr[i] += values[i];
  VecRestoreArray(rhs_, &rhs_ptr);
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
    MatSetValue(A_, rows[i], cols[i], vals[i], ADD_VALUES);
  MatAssemblyBegin(A_, MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FLUSH_ASSEMBLY);
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
                                  static_cast<int64_t>(sdm_.GetNumGhostDOFs(uk_man_)),
                                  ghids);
  }

  opensn::mpi_comm.barrier();
  log.Log() << "Done vector creation";
  opensn::mpi_comm.barrier();

  // Create KSP
  KSPCreate(opensn::mpi_comm, &ksp_);
  KSPSetOptionsPrefix(ksp_, name_.c_str());
  KSPSetType(ksp_, KSPCG);

  KSPSetTolerances(ksp_, 1.0e-50, options.residual_tolerance, 1.0e50, options.max_iters);

  // Set Pre-conditioner
  PC pc;
  KSPGetPC(ksp_, &pc);
  //  PCSetType(pc, PCGAMG);
  PCSetType(pc, PCHYPRE);

  PCHYPRESetType(pc, "boomeramg");
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
    PetscOptionsInsertString(nullptr, ("-" + name_ + option).c_str());

  PetscOptionsInsertString(nullptr, options.additional_options_string.c_str());

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp_);
}

void
DiffusionSolver::Solve(std::vector<double>& solution, bool use_initial_guess)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::Solve";
  Vec x;
  VecDuplicate(rhs_, &x);
  VecSet(x, 0.0);

  if (not use_initial_guess)
    KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
  else
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);

  KSPSetTolerances(
    ksp_, options.residual_tolerance, options.residual_tolerance, 1.0e50, options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr);

    double rhs_norm;
    VecNorm(rhs_, NORM_2, &rhs_norm);
    log.Log() << "RHS-norm " << rhs_norm;
  }

  if (use_initial_guess)
  {
    double* x_raw;
    VecGetArray(x, &x_raw);
    size_t k = 0;
    for (const auto& value : solution)
      x_raw[k++] = value;
    VecRestoreArray(x, &x_raw);
  }

  // Solve
  KSPSolve(ksp_, rhs_, x);

  // Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    log.Log() << "Solution-norm " << sol_norm;

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);

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
  VecDestroy(&x);
}

void
DiffusionSolver::Solve(Vec petsc_solution, bool use_initial_guess)
{
  const std::string fname = "acceleration::DiffusionMIPSolver::Solve";
  Vec x;
  VecDuplicate(rhs_, &x);
  VecSet(x, 0.0);

  if (not use_initial_guess)
    KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
  else
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);

  KSPSetTolerances(
    ksp_, options.residual_tolerance, options.residual_tolerance, 1.0e50, options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr);

    double rhs_norm;
    VecNorm(rhs_, NORM_2, &rhs_norm);
    log.Log() << "RHS-norm " << rhs_norm;
  }

  if (use_initial_guess)
  {
    VecCopy(petsc_solution, x);
  }

  // Solve
  KSPSolve(ksp_, rhs_, x);

  // Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    log.Log() << "Solution-norm " << sol_norm;

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);

    log.Log() << "Convergence Reason: " << GetPETScConvergedReasonstring(reason);
  }

  // Transfer petsc solution to vector
  VecCopy(x, petsc_solution);

  // Cleanup x
  VecDestroy(&x);
}

} // namespace opensn
