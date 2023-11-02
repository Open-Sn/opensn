#include "diffusion.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Initializes the diffusion solver. This involves creating the
 * sparse matrix with the appropriate sparsity pattern. Creating the
 * RHS vector. Creating the KSP solver. Setting the very specialized parameters
 * for Hypre's BooomerAMG. Note: `PCSetFromOptions` and
 * `KSPSetFromOptions` are called at the end. Therefore, any number of
 * additional PETSc options can be passed via the commandline.*/
void
lbs::acceleration::DiffusionSolver::Initialize()
{
  if (options.verbose) Chi::log.Log() << text_name_ << ": Initializing PETSc items";

  if (options.verbose)
    Chi::log.Log() << text_name_ << ": Global number of DOFs=" << num_global_dofs_;

  Chi::mpi.Barrier();
  Chi::log.Log() << "Sparsity pattern";
  Chi::mpi.Barrier();
  //============================================= Create Matrix
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm_.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, uk_man_);
  Chi::mpi.Barrier();
  Chi::log.Log() << "Done Sparsity pattern";
  Chi::mpi.Barrier();
  A_ = chi_math::PETScUtils::CreateSquareMatrix(num_local_dofs_, num_global_dofs_);
  chi_math::PETScUtils::InitMatrixSparsity(A_, nodal_nnz_in_diag, nodal_nnz_off_diag);
  Chi::mpi.Barrier();
  Chi::log.Log() << "Done matrix creation";
  Chi::mpi.Barrier();

  //============================================= Create RHS
  if (not requires_ghosts_)
    rhs_ = chi_math::PETScUtils::CreateVector(num_local_dofs_, num_global_dofs_);
  else
    rhs_ = chi_math::PETScUtils::CreateVectorWithGhosts(
      num_local_dofs_,
      num_global_dofs_,
      static_cast<int64_t>(sdm_.GetNumGhostDOFs(uk_man_)),
      sdm_.GetGhostDOFIndices(uk_man_));

  Chi::mpi.Barrier();
  Chi::log.Log() << "Done vector creation";
  Chi::mpi.Barrier();

  //============================================= Create KSP
  KSPCreate(PETSC_COMM_WORLD, &ksp_);
  KSPSetOptionsPrefix(ksp_, text_name_.c_str());
  KSPSetType(ksp_, KSPCG);

  KSPSetTolerances(ksp_, 1.e-50, options.residual_tolerance, 1.0e50, options.max_iters);

  //============================================= Set Pre-conditioner
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

  if (grid_.Attributes() & chi_mesh::DIMENSION_2)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.6");
  if (grid_.Attributes() & chi_mesh::DIMENSION_3)
    pc_options.emplace_back("pc_hypre_boomeramg_strong_threshold 0.8");

  for (const auto& option : pc_options)
    PetscOptionsInsertString(nullptr, ("-" + text_name_ + option).c_str());

  PetscOptionsInsertString(nullptr, options.additional_options_string.c_str());

  PCSetFromOptions(pc);
  KSPSetFromOptions(ksp_);
}