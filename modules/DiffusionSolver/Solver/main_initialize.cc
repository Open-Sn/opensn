#include "diffusion_solver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"

#include "chi_log.h"
#include "chi_mpi.h"

#include "utils/chi_timer.h"

// ###################################################################
/**Initializes the diffusion solver using the PETSc library.*/
int
chi_diffusion::Solver::Initialize(bool verbose)
{
  Chi::log.Log() << "\n"
                 << Chi::program_timer.GetTimeString() << " " << TextName()
                 << ": Initializing Diffusion solver ";
  this->verbose_info_ = verbose;

  if (not common_items_initialized_) InitializeCommonItems(); // Mostly boundaries

  chi::Timer t_init;
  t_init.Reset();

  auto sdm_string = basic_options_("discretization_method").StringValue();
  {
    using namespace chi_math::finite_element;
    if (sdm_string == "PWLC")
    {
      discretization_ =
        chi_math::spatial_discretization::PieceWiseLinearContinuous::New(*grid_ptr_);
      unknown_manager_.AddUnknown(chi_math::UnknownType::SCALAR);
    }
    else if (sdm_string == "PWLD_MIP")
    {
      discretization_ =
        chi_math::spatial_discretization::PieceWiseLinearDiscontinuous::New(*grid_ptr_);
      unknown_manager_.AddUnknown(chi_math::UnknownType::SCALAR);
    }
    else
      throw std::invalid_argument(TextName() + ": Invalid spatial discretization method, " +
                                  sdm_string + ", specified.");
  }

  unit_integrals_.clear();
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = discretization_->GetCellMapping(cell);
    unit_integrals_.insert(
      std::make_pair(cell.global_id_, UnitIntegralContainer::Make(cell_mapping)));
  }

  const auto ghost_global_ids = grid_ptr_->cells.GetGhostGlobalIDs();
  for (const uint64_t global_id : ghost_global_ids)
  {
    const auto& cell_mapping = discretization_->GetCellMapping(grid_ptr_->cells[global_id]);
    unit_integrals_.insert(std::make_pair(global_id, UnitIntegralContainer::Make(cell_mapping)));
  }

  MPI_Barrier(Chi::mpi.comm);
  auto& sdm = discretization_;

  //============================================= Get DOF counts
  local_dof_count_ = sdm->GetNumLocalDOFs(unknown_manager_);
  global_dof_count_ = sdm->GetNumGlobalDOFs(unknown_manager_);
  Chi::log.Log() << TextName() << ": Global number of DOFs=" << global_dof_count_;

  //================================================== Initialize discretization
  //                                                   method
  if (field_functions_.empty())
  {
    auto& sdm_ptr = discretization_;
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    using namespace chi_math;
    auto initial_field_function = std::make_shared<chi_physics::FieldFunctionGridBased>(
      text_name,                     // Text name
      sdm_ptr,                       // Spatial Discretization
      Unknown(UnknownType::SCALAR)); // Unknown/Variable

    field_functions_.push_back(initial_field_function);
    Chi::field_function_stack.push_back(initial_field_function);
  } // if not ff set

  //================================================== Determine nodal DOF
  Chi::log.Log() << "Building sparsity pattern.";
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm->BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, unknown_manager_);

  Chi::log.Log() << Chi::program_timer.GetTimeString() << " " << TextName()
                 << ": Diffusion Solver initialization time " << t_init.GetTime() / 1000.0
                 << std::endl;

  //================================================== Initialize x and b
  ierr_ = VecCreate(PETSC_COMM_WORLD, &x_);
  CHKERRQ(ierr_);
  ierr_ = PetscObjectSetName((PetscObject)x_, "Solution");
  CHKERRQ(ierr_);
  ierr_ = VecSetSizes(
    x_, static_cast<PetscInt>(local_dof_count_), static_cast<PetscInt>(global_dof_count_));
  CHKERRQ(ierr_);
  ierr_ = VecSetType(x_, VECMPI);
  CHKERRQ(ierr_);
  ierr_ = VecDuplicate(x_, &b_);
  CHKERRQ(ierr_);

  VecSet(x_, 0.0);
  VecSet(b_, 0.0);

  // ################################################## Create matrix
  ierr_ = MatCreate(PETSC_COMM_WORLD, &A_);
  CHKERRQ(ierr_);
  ierr_ = MatSetSizes(A_,
                      static_cast<PetscInt>(local_dof_count_),
                      static_cast<PetscInt>(local_dof_count_),
                      static_cast<PetscInt>(global_dof_count_),
                      static_cast<PetscInt>(global_dof_count_));
  CHKERRQ(ierr_);
  ierr_ = MatSetType(A_, MATMPIAIJ);
  CHKERRQ(ierr_);

  //================================================== Allocate matrix memory
  Chi::log.Log() << "Setting matrix preallocation.";
  MatMPIAIJSetPreallocation(A_, 0, nodal_nnz_in_diag.data(), 0, nodal_nnz_off_diag.data());
  MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetUp(A_);

  //================================================== Set up solver
  ierr_ = KSPCreate(PETSC_COMM_WORLD, &ksp_);
  ierr_ = KSPSetOperators(ksp_, A_, A_);
  ierr_ = KSPSetType(ksp_, KSPCG);

  //================================================== Set up preconditioner
  ierr_ = KSPGetPC(ksp_, &pc_);
  PCSetType(pc_, PCHYPRE);

  PCHYPRESetType(pc_, "boomeramg");

  //================================================== Setting Hypre parameters
  // The default HYPRE parameters used for polyhedra
  // seemed to have caused a lot of trouble for Slab
  // geometries. This section makes some custom options
  // per cell type
  auto first_cell = &grid_ptr_->local_cells[0];

  if (first_cell->Type() == chi_mesh::CellType::SLAB)
  {
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_P_max 4");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_grid_sweeps_coarse 1");

    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_interp_type ext+i");

    PetscOptionsInsertString(nullptr, "-options_left");
  }
  if (first_cell->Type() == chi_mesh::CellType::POLYGON)
  {
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_strong_threshold 0.6");

    // PetscOptionsInsertString(nullptr,"-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_P_max 4");

    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_interp_type ext+i");

    PetscOptionsInsertString(nullptr, "-options_left");
  }
  if (first_cell->Type() == chi_mesh::CellType::POLYHEDRON)
  {
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_strong_threshold 0.8");

    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_agg_nl 1");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_P_max 4");

    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_grid_sweeps_coarse 1");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_max_levels 25");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_relax_type_all symmetric-SOR/Jacobi");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_coarsen_type HMIS");
    PetscOptionsInsertString(nullptr, "-pc_hypre_boomeramg_interp_type ext+i");
  }
  PetscOptionsInsertString(nullptr, options_string_.c_str());
  PCSetFromOptions(pc_);

  //=================================== Set up monitor
  if (verbose) ierr_ = KSPMonitorSet(ksp_, &chi_diffusion::KSPMonitorAChiTech, nullptr, nullptr);

  KSPSetConvergenceTest(ksp_, &chi_diffusion::DiffusionConvergenceTestNPT, nullptr, nullptr);

  ierr_ = KSPSetTolerances(ksp_,
                           1.e-50,
                           basic_options_("residual_tolerance").FloatValue(),
                           1.0e50,
                           basic_options_("max_iters").IntegerValue());
  ierr_ = KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);

  return false;
}