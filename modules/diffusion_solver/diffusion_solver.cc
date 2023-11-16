#include "modules/diffusion_solver/diffusion_solver.h"
#include "framework/logging/log.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/cell_mappings/piecewise_linear_base_mapping.h"
#include "framework/physics/physics_material/physics_material.h"
#include "framework/physics/physics_material/material_property_scalar_value.h"
#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"
#include "framework/physics/field_function/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/mpi/mpi.h"
#include "framework/utils/timer.h"

namespace opensn
{
namespace diffusion
{

Solver::Solver(const std::string& in_solver_name)
  : opensn::Solver(in_solver_name,
                   {{"discretization_method", std::string("None")},
                    {"max_iters", int64_t(500)},
                    {"residual_tolerance", 1.0e-8},
                    {"property_map_D", int64_t(0)},
                    {"property_map_q", int64_t(1)},
                    {"property_map_sigma", int64_t(2)}})
{
}

Solver::~Solver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
  KSPDestroy(&ksp_);
}

void
Solver::GetMaterialProperties(const Cell& cell,
                              int cell_dofs,
                              std::vector<double>& diffCoeff,
                              std::vector<double>& sourceQ,
                              std::vector<double>& sigmaa,
                              int group,
                              int moment)
{
  uint64_t cell_glob_index = cell.global_id_;
  bool cell_is_local = (cell.partition_id_ == opensn::mpi.location_id);
  uint64_t cell_local_id = cell.local_id_;
  int mat_id = cell.material_id_;

  if (mat_id < 0)
  {
    log.Log0Error() << "Cell encountered with no material id. ";
    Exit(EXIT_FAILURE);
  }

  if (mat_id >= Chi::material_stack.size())
  {
    log.Log0Error() << "Cell encountered with material id pointing to "
                       "non-existing material.";
    Exit(EXIT_FAILURE);
  }

  auto property_map_D = basic_options_("property_map_D").IntegerValue();
  auto property_map_q = basic_options_("property_map_q").IntegerValue();
  auto property_map_sigma = basic_options_("property_map_sigma").IntegerValue();

  auto material = Chi::GetStackItemPtr(Chi::material_stack, mat_id, __FUNCTION__);

  // Process material properties
  diffCoeff.resize(cell_dofs, 1.0);
  sourceQ.resize(cell_dofs, 0.0);
  sigmaa.resize(cell_dofs, 0.0);

  // REGULAR MATERIAL
  if (material_mode_ == DIFFUSION_MATERIALS_REGULAR)
  {
    // We absolutely need the diffusion coefficient so process error
    if ((property_map_D < 0) || (property_map_D >= material->properties_.size()))
    {
      log.Log0Error() << "Solver diffusion coefficient mapped to property index " << property_map_D
                      << " is not a valid index for material \"" << material->name_ << "\" id "
                      << mat_id;
      Exit(EXIT_FAILURE);
    }

    // For now, we can only support scalar values so lets check that
    if (std::dynamic_pointer_cast<ScalarValue>(material->properties_[property_map_D]))
    {
      diffCoeff.assign(cell_dofs, material->properties_[property_map_D]->GetScalarValue());
    }
    else
    {
      log.Log0Error() << "Solver diffusion coefficient mapped to property index " << property_map_D
                      << " is not a valid property type"
                      << " for material \"" << material->name_ << "\" id " << mat_id
                      << ". Currently SCALAR_VALUE and THERMAL_CONDUCTIVITY are the "
                      << "only supported types.";
      Exit(EXIT_FAILURE);
    }

    if ((property_map_q < material->properties_.size()) && (property_map_q >= 0))
    {
      if (std::dynamic_pointer_cast<ScalarValue>(material->properties_[property_map_q]))
      {
        sourceQ.assign(cell_dofs, material->properties_[property_map_q]->GetScalarValue());
      }
      else
      {
        log.Log0Error() << "Source value mapped to property index " << property_map_q
                        << " is not a valid property type"
                        << " for material \"" << material->name_ << "\" id " << mat_id
                        << ". Currently SCALAR_VALUE is the "
                        << "only supported type.";
        Exit(EXIT_FAILURE);
      }
    }

    if (not((property_map_sigma < 0) || (property_map_sigma >= material->properties_.size())))
    {
      sigmaa.assign(cell_dofs, material->properties_[property_map_sigma]->GetScalarValue());
    }
  } // regular

  // TRANSPORT XS D
  // TRANSPORT XS SIGA
  // SCALAR       Q
  else if (material_mode_ == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR)
  {
    // Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p = 0; p < material->properties_.size(); p++)
    {
      if (std::dynamic_pointer_cast<MultiGroupXS>(material->properties_[p]))
      {
        auto xs = std::static_pointer_cast<MultiGroupXS>(material->properties_[p]);

        diffCoeff.assign(cell_dofs, xs->DiffusionCoefficient()[group]);
        sigmaa.assign(cell_dofs, xs->SigmaRemoval()[group]);
        transportxs_found = true;
      }
    } // for properties

    if (!transportxs_found)
    {
      log.LogAllError() << "Diffusion Solver: Material encountered with no tranport xs"
                           " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      Exit(EXIT_FAILURE);
    }

    // Setting Q
    if ((property_map_q < material->properties_.size()) && (property_map_q >= 0))
    {
      if (std::dynamic_pointer_cast<ScalarValue>(material->properties_[property_map_q]))
      {
        sourceQ.assign(cell_dofs, material->properties_[property_map_q]->GetScalarValue());
      }
      else
      {
        log.Log0Error() << "Source value mapped to property index " << property_map_q
                        << " is not a valid property type"
                        << " for material \"" << material->name_ << "\" id " << mat_id
                        << ". Currently SCALAR_VALUE is the "
                        << "only supported type.";
        Exit(EXIT_FAILURE);
      }
    }
  } // transport xs TTR
  else
  {
    log.Log0Error() << "Diffusion Solver: Invalid material mode.";
    Exit(EXIT_FAILURE);
  }
}

void
Solver::UpdateFieldFunctions()
{
  log.LogAll() << "Updating field functions" << std::endl;
  auto& ff = *field_functions_.front();
  const auto& OneDofPerNode = discretization_->UNITARY_UNKNOWN_MANAGER;

  std::vector<double> data_vector;
  discretization_->LocalizePETScVector(x_, data_vector, OneDofPerNode);

  ff.UpdateFieldVector(data_vector);
}

void
Solver::InitializeCommonItems()
{
  const std::string fname = "Solver::InitializeCommonItems";
  grid_ptr_ = GetCurrentHandler().GetGrid();

  if (grid_ptr_ == nullptr) throw std::logic_error(fname + " No grid defined.");

  auto globl_unique_bndry_ids = grid_ptr_->GetDomainUniqueBoundaryIDs();

  const auto& grid_boundary_id_map = grid_ptr_->GetBoundaryIDMap();
  for (uint64_t bndry_id : globl_unique_bndry_ids)
  {
    if (grid_boundary_id_map.count(bndry_id) == 0)
      throw std::logic_error(fname + ": Boundary id " + std::to_string(bndry_id) +
                             " does not have a name-assignment.");

    const auto& bndry_name = grid_boundary_id_map.at(bndry_id);
    if (boundary_preferences_.find(bndry_name) != boundary_preferences_.end())
    {
      BoundaryInfo bndry_info = boundary_preferences_.at(bndry_name);
      auto& bndry_vals = bndry_info.second;
      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting:
        {
          boundaries_.insert(std::make_pair(bndry_id, new BoundaryReflecting()));
          log.Log() << "Boundary \"" << bndry_name << "\" set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1, 0.0);
          boundaries_.insert(std::make_pair(bndry_id, new BoundaryDirichlet(bndry_vals[0])));
          log.Log() << "Boundary \"" << bndry_name << "\" set to dirichlet.";
          break;
        }
        case BoundaryType::Robin:
        {
          if (bndry_vals.size() < 3) bndry_vals.resize(3, 0.0);
          boundaries_.insert(std::make_pair(
            bndry_id, new BoundaryRobin(bndry_vals[0], bndry_vals[1], bndry_vals[2])));
          log.Log() << "Boundary \"" << bndry_name << "\" set to robin.";
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries_.insert(std::make_pair(bndry_id, new BoundaryRobin(0.25, 0.5, 0.0)));
          log.Log() << "Boundary \"" << bndry_name << "\" set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size() < 3) bndry_vals.resize(3, 0.0);
          boundaries_.insert(std::make_pair(
            bndry_id, new BoundaryRobin(bndry_vals[0], bndry_vals[1], bndry_vals[2])));
          log.Log() << "Boundary \"" << bndry_name << "\" set to neumann.";
          break;
        }
      } // switch boundary type
    }
    else
    {
      boundaries_.insert(std::make_pair(bndry_id, new BoundaryDirichlet()));
      log.Log0Verbose1() << "No boundary preference found for boundary \"" << bndry_name
                         << "\" Dirichlet boundary added with zero boundary value.";
    }
  } // for neighbor_id_

  common_items_initialized_ = true;
}

int
Solver::Initialize(bool verbose)
{
  log.Log() << "\n"
            << program_timer.GetTimeString() << " " << TextName()
            << ": Initializing Diffusion solver ";
  this->verbose_info_ = verbose;

  if (not common_items_initialized_) InitializeCommonItems(); // Mostly boundaries

  Timer t_init;
  t_init.Reset();

  auto sdm_string = basic_options_("discretization_method").StringValue();
  {
    if (sdm_string == "PWLC")
    {
      discretization_ = PieceWiseLinearContinuous::New(*grid_ptr_);
      unknown_manager_.AddUnknown(UnknownType::SCALAR);
    }
    else if (sdm_string == "PWLD_MIP")
    {
      discretization_ = PieceWiseLinearDiscontinuous::New(*grid_ptr_);
      unknown_manager_.AddUnknown(UnknownType::SCALAR);
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

  MPI_Barrier(mpi.comm);
  auto& sdm = discretization_;

  // Get DOF counts
  local_dof_count_ = sdm->GetNumLocalDOFs(unknown_manager_);
  global_dof_count_ = sdm->GetNumGlobalDOFs(unknown_manager_);
  log.Log() << TextName() << ": Global number of DOFs=" << global_dof_count_;

  // Initialize discretization method
  if (field_functions_.empty())
  {
    auto& sdm_ptr = discretization_;
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    auto initial_field_function =
      std::make_shared<FieldFunctionGridBased>(text_name, sdm_ptr, Unknown(UnknownType::SCALAR));

    field_functions_.push_back(initial_field_function);
    Chi::field_function_stack.push_back(initial_field_function);
  } // if not ff set

  // Determine nodal DOF
  log.Log() << "Building sparsity pattern.";
  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm->BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, unknown_manager_);

  log.Log() << program_timer.GetTimeString() << " " << TextName()
            << ": Diffusion Solver initialization time " << t_init.GetTime() / 1000.0 << std::endl;

  // Initialize x and b
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

  // Create matrix
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

  // Allocate matrix memory
  log.Log() << "Setting matrix preallocation.";
  MatMPIAIJSetPreallocation(A_, 0, nodal_nnz_in_diag.data(), 0, nodal_nnz_off_diag.data());
  MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetUp(A_);

  // Set up solver
  ierr_ = KSPCreate(PETSC_COMM_WORLD, &ksp_);
  ierr_ = KSPSetOperators(ksp_, A_, A_);
  ierr_ = KSPSetType(ksp_, KSPCG);

  // Set up preconditioner
  ierr_ = KSPGetPC(ksp_, &pc_);
  PCSetType(pc_, PCHYPRE);

  PCHYPRESetType(pc_, "boomeramg");

  // Setting Hypre parameters
  // The default HYPRE parameters used for polyhedra
  // seemed to have caused a lot of trouble for Slab
  // geometries. This section makes some custom options
  // per cell type
  auto first_cell = &grid_ptr_->local_cells[0];

  if (first_cell->Type() == CellType::SLAB)
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
  if (first_cell->Type() == CellType::POLYGON)
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
  if (first_cell->Type() == CellType::POLYHEDRON)
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

  // Set up monitor
  if (verbose) ierr_ = KSPMonitorSet(ksp_, &KSPMonitorAChiTech, nullptr, nullptr);

  KSPSetConvergenceTest(ksp_, &DiffusionConvergenceTestNPT, nullptr, nullptr);

  ierr_ = KSPSetTolerances(ksp_,
                           1.e-50,
                           basic_options_("residual_tolerance").FloatValue(),
                           1.0e50,
                           basic_options_("max_iters").IntegerValue());
  ierr_ = KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);

  return false;
}

int
Solver::ExecuteS(bool suppress_assembly, bool suppress_solve)
{
  t_assembly_.Reset();

  if (Chi::material_stack.empty())
  {
    log.Log0Error() << "No materials added to simulation. Add materials.";
    exit(0);
  }

  VecSet(x_, 0.0);
  VecSet(b_, 0.0);

  if (!suppress_assembly)
    log.Log() << program_timer.GetTimeString() << " " << TextName() << ": Assembling A locally";

  // Loop over locally owned cells
  auto fem_method = basic_options_("discretization_method").StringValue();
  if (fem_method == "PWLC")
  {
    if (!suppress_assembly)
      for (auto& cell : grid_ptr_->local_cells)
        CFEM_Assemble_A_and_b(cell, gi_);
    else {}
  }
  else if (fem_method == "PWLD_MIP")
  {
    if (!suppress_assembly)
      for (auto& cell : grid_ptr_->local_cells)
        PWLD_Assemble_A_and_b(cell, gi_);
    else
      for (auto& cell : grid_ptr_->local_cells)
        PWLD_Assemble_b(cell, gi_);
  }
  else
  {
    log.Log() << "Diffusion Solver: Finite Element Discretization "
                 "method not specified.";
    Exit(EXIT_FAILURE);
  }

  if (!suppress_assembly)
    log.Log() << program_timer.GetTimeString() << " " << TextName()
              << ": Done Assembling A locally";
  MPI_Barrier(mpi.comm);

  // Call matrix assembly
  if (verbose_info_ || log.GetVerbosity() >= Logger::LOG_0VERBOSE_1)
    log.Log() << program_timer.GetTimeString() << " " << TextName()
              << ": Communicating matrix assembly";

  if (!suppress_assembly)
  {
    log.Log() << program_timer.GetTimeString() << " " << TextName() << ": Assembling A globally";
    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);

    // Matrix symmetry check
    //    PetscBool is_symmetric;
    //    ierr = MatIsSymmetric(A,1.0e-4,&is_symmetric);
    //    if (!is_symmetric)
    //    {
    //      chi::log.Log0Warning()
    //        << "Assembled matrix is not symmetric";
    //    }

    // Matrix diagonal check
    log.Log() << program_timer.GetTimeString() << " " << TextName() << ": Diagonal check";
    PetscBool missing_diagonal;
    PetscInt row;
    MatMissingDiagonal(A_, &missing_diagonal, &row);
    if (missing_diagonal)
      log.LogAllError() << program_timer.GetTimeString() << " " << TextName()
                        << ": Missing diagonal detected";

    // Matrix sparsity info
    MatInfo info;
    ierr_ = MatGetInfo(A_, MAT_GLOBAL_SUM, &info);

    log.Log() << "Number of mallocs used = " << info.mallocs
              << "\nNumber of non-zeros allocated = " << info.nz_allocated
              << "\nNumber of non-zeros used = " << info.nz_used
              << "\nNumber of unneeded non-zeros = " << info.nz_unneeded;
  }
  if (verbose_info_ || log.GetVerbosity() >= Logger::LOG_0VERBOSE_1)
    log.Log() << program_timer.GetTimeString() << " " << TextName() << ": Assembling x and b";
  VecAssemblyBegin(x_);
  VecAssemblyEnd(x_);
  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

  time_assembly_ = t_assembly_.GetTime() / 1000.0;

  // Execute solve
  if (suppress_solve)
  {
    log.Log() << program_timer.GetTimeString() << " " << TextName()
              << ": Setting up solver and preconditioner\n";
    PCSetUp(pc_);
    KSPSetUp(ksp_);
  }
  else
  {
    if (verbose_info_ || log.GetVerbosity() >= Logger::LOG_0VERBOSE_1)
      log.Log() << program_timer.GetTimeString() << " " << TextName() << ": Solving system\n";
    t_solve_.Reset();
    PCSetUp(pc_);
    KSPSetUp(ksp_);
    KSPSolve(ksp_, b_, x_);
    time_solve_ = t_solve_.GetTime() / 1000.0;

    // Populate field vector
    //    if (fem_method == "PWLD_MIP" or fem_method == "PWLD_MIP_GAGG")
    //    {
    //      const double* x_ref;
    //      VecGetArrayRead(x,&x_ref);
    //
    //      for (int i=0; i < local_dof_count; i++)
    //        pwld_phi_local[i] = x_ref[i];
    //
    //      VecRestoreArrayRead(x,&x_ref);
    //    }

    // Get convergence reason
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);
    if (verbose_info_ || reason != KSP_CONVERGED_RTOL)
      log.Log() << "Convergence reason: " << GetPETScConvergedReasonstring(reason);

    // Location wise view
    if (opensn::mpi.location_id == 0)
    {
      int64_t its;
      ierr_ = KSPGetIterationNumber(ksp_, &its);
      log.Log() << program_timer.GetTimeString() << " " << TextName() << "[g=" << gi_ << "-"
                << gi_ + G_ - 1 << "]: Number of iterations =" << its;

      if (verbose_info_ || log.GetVerbosity() >= Logger::LOG_0VERBOSE_1)
      {
        log.Log() << "Timing:";
        log.Log() << "Assembling the matrix: " << time_assembly_;
        log.Log() << "Solving the system   : " << time_solve_;
      }
    }

    UpdateFieldFunctions();

    if (verbose_info_ || log.GetVerbosity() >= Logger::LOG_0VERBOSE_1)
      log.Log() << "Diffusion Solver execution completed!\n";
  } // if not suppressed solve

  return 0;
}

void
Solver::CFEM_Assemble_A_and_b(Cell& cell, int group)
{
  auto pwl_sdm = std::static_pointer_cast<PieceWiseLinearContinuous>(this->discretization_);
  const auto& fe_intgrl_values = unit_integrals_.at(cell.global_id_);

  size_t num_nodes = fe_intgrl_values.NumNodes();

  // Process material id
  std::vector<double> D(num_nodes, 1.0);
  std::vector<double> q(num_nodes, 1.0);
  std::vector<double> siga(num_nodes, 0.0);

  GetMaterialProperties(cell, num_nodes, D, q, siga, group);

  // Init cell matrix info
  typedef std::vector<double> Row;
  typedef std::vector<Row> Matrix;

  Matrix cell_matrix;
  std::vector<double> cell_rhs;

  cell_matrix.resize(num_nodes, Row(num_nodes, 0.0));
  cell_rhs.resize(num_nodes, 0.0);

  std::vector<int64_t> dof_global_row_ind(num_nodes, -1);
  std::vector<int64_t> dof_global_col_ind(num_nodes, -1);

  // Loop over DOFs
  for (int i = 0; i < num_nodes; i++)
  {
    dof_global_row_ind[i] = pwl_sdm->MapDOF(cell, i);

    for (int j = 0; j < num_nodes; j++)
    {
      double mat_entry = D[j] * fe_intgrl_values.IntV_gradShapeI_gradShapeJ(i, j) +
                         siga[j] * fe_intgrl_values.IntV_shapeI_shapeJ(i, j);

      cell_matrix[i][j] = mat_entry;
    } // for j

    // Develop RHS entry
    cell_rhs[i] = q[i] * fe_intgrl_values.IntV_shapeI(i);
  } // for i
  dof_global_col_ind = dof_global_row_ind;

  // Apply Dirichlet,Vacuum, Neumann and Robin BCs
  // Dirichlets are just collected
  std::vector<int> dirichlet_count(num_nodes, 0);
  std::vector<double> dirichlet_value(num_nodes, 0.0);
  for (int f = 0; f < cell.faces_.size(); f++)
  {
    if (not cell.faces_[f].has_neighbor_)
    {
      uint64_t ir_boundary_index = cell.faces_[f].neighbor_id_;
      auto ir_boundary_type = boundaries_.at(ir_boundary_index)->type_;

      if (ir_boundary_type == BoundaryType::Dirichlet)
      {
        auto& dirichlet_bndry = (BoundaryDirichlet&)*boundaries_.at(ir_boundary_index);

        int num_face_dofs = cell.faces_[f].vertex_ids_.size();
        for (int fi = 0; fi < num_face_dofs; fi++)
        {
          int i = fe_intgrl_values.FaceDofMapping(f, fi);
          dirichlet_count[i] += 1;
          dirichlet_value[i] += dirichlet_bndry.boundary_value;
        }
      }

      if (ir_boundary_type == BoundaryType::Robin)
      {
        auto& robin_bndry = (BoundaryRobin&)*boundaries_.at(ir_boundary_index);

        std::cout << robin_bndry.a << " " << robin_bndry.b << " " << robin_bndry.f << std::endl;

        int num_face_dofs = cell.faces_[f].vertex_ids_.size();
        for (int fi = 0; fi < num_face_dofs; fi++)
        {
          int i = fe_intgrl_values.FaceDofMapping(f, fi);

          for (int fj = 0; fj < num_face_dofs; fj++)
          {
            int j = fe_intgrl_values.FaceDofMapping(f, fj);

            double aij = robin_bndry.a * fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);
            aij /= robin_bndry.b;

            cell_matrix[i][j] += aij;
          } // for fj

          double aii = robin_bndry.f * fe_intgrl_values.IntS_shapeI(f, i);
          aii /= robin_bndry.b;

          cell_matrix[i][i] += aii;
        } // for fi
      }   // if robin

    } // if boundary

  } // for face

  // Apply dirichlet BCs
  // Compute average dirichlet value
  for (int i = 0; i < num_nodes; ++i)
    dirichlet_value[i] /= (dirichlet_count[i] > 0) ? dirichlet_count[i] : 1;

  for (int i = 0; i < num_nodes; ++i)
  {
    if (dirichlet_count[i] > 0)
    {
      cell_matrix[i].clear();
      cell_matrix[i] = std::vector<double>(num_nodes, 0.0);
      cell_matrix[i][i] = 1.0;
      int ir = dof_global_col_ind[i];
      MatSetValue(A_, ir, ir, 1.0, ADD_VALUES);
      dof_global_col_ind[i] = -1;
      cell_rhs[i] = dirichlet_value[i];
    }
    else
    {
      for (int j = 0; j < num_nodes; ++j)
      {
        if (dirichlet_count[j] > 0)
        {
          cell_rhs[i] -= cell_matrix[i][j] * dirichlet_value[j];
          cell_matrix[i][j] = 0.0;
        }
      }
    }
  }

  // Make contiguous copy of matrix
  std::vector<double> cell_matrix_cont(num_nodes * num_nodes, 0.0);
  int n = 0;
  for (int i = 0; i < num_nodes; ++i)
    for (int j = 0; j < num_nodes; ++j)
      cell_matrix_cont[n++] = cell_matrix[i][j];

  // Add to global
  MatSetValues(A_,
               num_nodes,
               dof_global_row_ind.data(),
               num_nodes,
               dof_global_col_ind.data(),
               cell_matrix_cont.data(),
               ADD_VALUES);

  VecSetValues(b_, num_nodes, dof_global_row_ind.data(), cell_rhs.data(), ADD_VALUES);

  VecSetValues(x_, num_nodes, dof_global_row_ind.data(), dirichlet_value.data(), INSERT_VALUES);
}

void
Solver::PWLD_Assemble_A_and_b(const Cell& cell, int component)
{
  auto pwl_sdm = std::static_pointer_cast<PieceWiseLinearDiscontinuous>(this->discretization_);
  const auto& fe_intgrl_values = unit_integrals_.at(cell.global_id_);

  size_t num_nodes = fe_intgrl_values.NumNodes();

  // Process material id
  int mat_id = cell.material_id_;

  std::vector<double> D(num_nodes, 1.0);
  std::vector<double> q(num_nodes, 1.0);
  std::vector<double> siga(num_nodes, 0.0);

  GetMaterialProperties(cell, num_nodes, D, q, siga, component);

  // Loop over DOFs
  for (int i = 0; i < num_nodes; i++)
  {
    int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);
    double rhsvalue = 0.0;

    // Develop matrix entry
    for (int j = 0; j < num_nodes; j++)
    {
      int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);

      double jr_mat_entry = D[j] * fe_intgrl_values.IntV_gradShapeI_gradShapeJ(i, j);

      jr_mat_entry += siga[j] * fe_intgrl_values.IntV_shapeI_shapeJ(i, j);

      MatSetValue(A_, ir, jr, jr_mat_entry, ADD_VALUES);

      rhsvalue += q[j] * fe_intgrl_values.IntV_shapeI_shapeJ(i, j);
    } // for j

    // Apply RHS entry
    VecSetValue(b_, ir, rhsvalue, ADD_VALUES);

  } // for i

  // Loop over faces
  int num_faces = cell.faces_.size();
  for (unsigned int f = 0; f < num_faces; f++)
  {
    auto& face = cell.faces_[f];

    // Get face normal
    Vector3 n = face.normal_;

    int num_face_dofs = face.vertex_ids_.size();

    if (face.has_neighbor_)
    {
      const auto& adj_cell = grid_ptr_->cells[face.neighbor_id_];
      const auto& adj_fe_intgrl_values = unit_integrals_.at(adj_cell.global_id_);

      // Get the current map to the adj cell's face
      unsigned int fmap = MapCellFace(cell, adj_cell, f);

      // Compute penalty coefficient
      double hp = HPerpendicular(adj_cell, adj_fe_intgrl_values, fmap);
      double hm = HPerpendicular(cell, fe_intgrl_values, f);

      std::vector<double> adj_D, adj_Q, adj_sigma;

      GetMaterialProperties(
        adj_cell, adj_fe_intgrl_values.NumNodes(), adj_D, adj_Q, adj_sigma, component);

      // Compute surface average D
      double D_avg = 0.0;
      double intS = 0.0;
      for (int fi = 0; fi < num_face_dofs; fi++)
      {
        int i = fe_intgrl_values.FaceDofMapping(f, fi);
        D_avg += D[i] * fe_intgrl_values.IntS_shapeI(f, i);
        intS += fe_intgrl_values.IntS_shapeI(f, i);
      }
      D_avg /= intS;

      // Compute surface average D_adj
      double adj_D_avg = 0.0;
      double adj_intS = 0.0;
      for (int fi = 0; fi < num_face_dofs; fi++)
      {
        int i = fe_intgrl_values.FaceDofMapping(f, fi);
        int imap = MapCellLocalNodeIDFromGlobalID(adj_cell, cell.vertex_ids_[i]);
        adj_D_avg += adj_D[imap] * adj_fe_intgrl_values.IntS_shapeI(fmap, imap);
        adj_intS += adj_fe_intgrl_values.IntS_shapeI(fmap, imap);
      }
      adj_D_avg /= adj_intS;

      // Compute kappa
      double kappa = 1.0;
      if (cell.Type() == CellType::SLAB) kappa = fmax(2.0 * (adj_D_avg / hp + D_avg / hm), 0.25);
      if (cell.Type() == CellType::POLYGON) kappa = fmax(2.0 * (adj_D_avg / hp + D_avg / hm), 0.25);
      if (cell.Type() == CellType::POLYHEDRON)
        kappa = fmax(4.0 * (adj_D_avg / hp + D_avg / hm), 0.25);

      // Assembly penalty terms
      for (int fi = 0; fi < num_face_dofs; fi++)
      {
        int i = fe_intgrl_values.FaceDofMapping(f, fi);
        int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);

        for (int fj = 0; fj < num_face_dofs; fj++)
        {
          int j = fe_intgrl_values.FaceDofMapping(f, fj);
          int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);
          int jmap = MapCellLocalNodeIDFromGlobalID(adj_cell, face.vertex_ids_[fj]);
          int jrmap = pwl_sdm->MapDOF(adj_cell, jmap, unknown_manager_, 0, component);

          double aij = kappa * fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);

          MatSetValue(A_, ir, jr, aij, ADD_VALUES);
          MatSetValue(A_, ir, jrmap, -aij, ADD_VALUES);
        } // for fj

      } // for fi

      // Assemble gradient terms
      // For the following comments we use the notation:
      // Dk = 0.5* n dot nabla bk

      // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
      for (int i = 0; i < fe_intgrl_values.NumNodes(); i++)
      {
        int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);

        for (int fj = 0; fj < num_face_dofs; fj++)
        {
          int j = fe_intgrl_values.FaceDofMapping(f, fj);
          int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);
          int jmap = MapCellLocalNodeIDFromGlobalID(adj_cell, face.vertex_ids_[fj]);
          int jrmap = pwl_sdm->MapDOF(adj_cell, jmap, unknown_manager_, 0, component);

          double aij = -0.5 * D_avg * n.Dot(fe_intgrl_values.IntS_shapeI_gradshapeJ(f, j, i));

          MatSetValue(A_, ir, jr, aij, ADD_VALUES);
          MatSetValue(A_, ir, jrmap, -aij, ADD_VALUES);
        } // for fj
      }   // for i

      // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
      for (int fi = 0; fi < num_face_dofs; fi++)
      {
        int i = fe_intgrl_values.FaceDofMapping(f, fi);
        int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);
        int imap = MapCellLocalNodeIDFromGlobalID(adj_cell, face.vertex_ids_[fi]);
        int irmap = pwl_sdm->MapDOF(adj_cell, imap, unknown_manager_, 0, component);

        for (int j = 0; j < fe_intgrl_values.NumNodes(); j++)
        {
          int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);

          double aij = -0.5 * D_avg * n.Dot(fe_intgrl_values.IntS_shapeI_gradshapeJ(f, i, j));

          MatSetValue(A_, ir, jr, aij, ADD_VALUES);
          MatSetValue(A_, irmap, jr, -aij, ADD_VALUES);
        } // for j
      }   // for fi

    } // if not bndry
    else
    {
      uint64_t ir_boundary_index = face.neighbor_id_;
      auto ir_boundary_type = boundaries_.at(ir_boundary_index)->type_;

      if (ir_boundary_type == BoundaryType::Dirichlet)
      {
        auto dc_boundary = (BoundaryDirichlet&)*boundaries_.at(ir_boundary_index);

        // Compute penalty coefficient
        double hm = HPerpendicular(cell, fe_intgrl_values, f);

        // Compute surface average D
        double D_avg = 0.0;
        double intS = 0.0;
        for (int fi = 0; fi < num_face_dofs; fi++)
        {
          int i = fe_intgrl_values.FaceDofMapping(f, fi);
          D_avg += D[i] * fe_intgrl_values.IntS_shapeI(f, i);
          intS += fe_intgrl_values.IntS_shapeI(f, i);
        }
        D_avg /= intS;

        double kappa = 1.0;
        if (cell.Type() == CellType::SLAB) kappa = fmax(4.0 * (D_avg / hm), 0.25);
        if (cell.Type() == CellType::POLYGON) kappa = fmax(4.0 * (D_avg / hm), 0.25);
        if (cell.Type() == CellType::POLYHEDRON) kappa = fmax(8.0 * (D_avg / hm), 0.25);

        // Assembly penalty terms
        for (int fi = 0; fi < num_face_dofs; fi++)
        {
          int i = fe_intgrl_values.FaceDofMapping(f, fi);
          int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);

          for (int fj = 0; fj < num_face_dofs; fj++)
          {
            int j = fe_intgrl_values.FaceDofMapping(f, fj);
            int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);

            double aij = kappa * fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);

            MatSetValue(A_, ir, jr, aij, ADD_VALUES);
            VecSetValue(b_, ir, aij * dc_boundary.boundary_value, ADD_VALUES);
          } // for fj
        }   // for fi

        // -Di^- bj^- and
        // -Dj^- bi^-
        for (int i = 0; i < fe_intgrl_values.NumNodes(); i++)
        {
          int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);

          for (int j = 0; j < fe_intgrl_values.NumNodes(); j++)
          {
            int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);

            double gij = n.Dot(fe_intgrl_values.IntS_shapeI_gradshapeJ(f, i, j) +
                               fe_intgrl_values.IntS_shapeI_gradshapeJ(f, j, i));
            double aij = -0.5 * D_avg * gij;

            MatSetValue(A_, ir, jr, aij, ADD_VALUES);
            VecSetValue(b_, ir, aij * dc_boundary.boundary_value, ADD_VALUES);
          } // for j
        }   // for i
      }     // Dirichlet
      else if (ir_boundary_type == BoundaryType::Robin)
      {
        auto robin_bndry = (BoundaryRobin&)*boundaries_.at(ir_boundary_index);

        for (int fi = 0; fi < num_face_dofs; fi++)
        {
          int i = fe_intgrl_values.FaceDofMapping(f, fi);
          int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);

          for (int fj = 0; fj < num_face_dofs; fj++)
          {
            int j = fe_intgrl_values.FaceDofMapping(f, fj);
            int jr = pwl_sdm->MapDOF(cell, j, unknown_manager_, 0, component);

            double aij = robin_bndry.a * fe_intgrl_values.IntS_shapeI_shapeJ(f, i, j);
            aij /= robin_bndry.b;

            MatSetValue(A_, ir, jr, aij, ADD_VALUES);
          } // for fj

          double bi = robin_bndry.f * fe_intgrl_values.IntS_shapeI(f, i);
          bi /= robin_bndry.b;

          VecSetValue(b_, ir, bi, ADD_VALUES);
        } // for fi
      }   // robin
    }
  } // for f
}

void
Solver::PWLD_Assemble_b(const Cell& cell, int component)
{
  auto pwl_sdm = std::static_pointer_cast<PieceWiseLinearDiscontinuous>(this->discretization_);
  const auto& fe_intgrl_values = unit_integrals_.at(cell.global_id_);

  size_t num_nodes = fe_intgrl_values.NumNodes();

  // Process material id
  int mat_id = cell.material_id_;

  std::vector<double> D(num_nodes, 1.0);
  std::vector<double> q(num_nodes, 1.0);
  std::vector<double> siga(num_nodes, 1.0);

  GetMaterialProperties(cell, num_nodes, D, q, siga, component);

  // Loop over DOFs
  for (int i = 0; i < num_nodes; i++)
  {
    int ir = pwl_sdm->MapDOF(cell, i, unknown_manager_, 0, component);

    // Develop rhs entry
    double rhsvalue = 0.0;
    for (int j = 0; j < num_nodes; j++)
      rhsvalue += q[j] * fe_intgrl_values.IntV_shapeI_shapeJ(i, j);

    // Apply RHS entry
    VecSetValue(b_, ir, rhsvalue, ADD_VALUES);

  } // for i
}

double
Solver::HPerpendicular(const Cell& cell,
                       const UnitIntegralContainer& fe_intgrl_values,
                       unsigned int f)
{
  double hp = 1.0;

  size_t Nf = cell.faces_.size();
  size_t Nv = cell.vertex_ids_.size();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == CellType::SLAB)
  {
    const auto& v0 = grid_ptr_->vertices[cell.vertex_ids_[0]];
    const auto& v1 = grid_ptr_->vertices[cell.vertex_ids_[1]];

    hp = (v1 - v0).Norm() / 2.0;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == CellType::POLYGON)
  {
    //    Nv = 4;
    const CellFace& face = cell.faces_[f];

    uint64_t v0i = face.vertex_ids_[0];
    uint64_t v1i = face.vertex_ids_[1];

    const auto& v0 = grid_ptr_->vertices[v0i];
    const auto& v1 = grid_ptr_->vertices[v1i];

    double perimeter = (v1 - v0).Norm();

    double area = 0.0;
    for (int i = 0; i < fe_intgrl_values.NumNodes(); i++)
      area += fe_intgrl_values.IntV_shapeI(i);

    hp = area / perimeter;

    //    if (Nv == 3)
    //      hp = 2*area/perimeter;
    //    else if (Nv == 4)
    //      hp = area/perimeter;
    //    else //Nv > 4
    //    {
    //      if (Nv%2 == 0)
    //        hp = 4*area/perimeter;
    //      else
    //      {
    //        hp = 2*area/perimeter;
    //        hp += sqrt(2*area / Nv*sin(2*M_PI/Nv));
    //      }
    //    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == CellType::POLYHEDRON)
  {
    double volume = 0.0;
    for (int i = 0; i < fe_intgrl_values.NumNodes(); i++)
      volume += fe_intgrl_values.IntV_shapeI(i);

    double area = 0.0;
    for (int fr = 0; fr < Nf; fr++)
      for (int i = 0; i < Nv; i++)
        area += fe_intgrl_values.IntS_shapeI(fr, i);

    if (Nf == 4) // Tet
      hp = 3 * volume / area;
    else if (Nf == 6 && Nv == 8) // Hex
      hp = volume / area;
    else // Polyhedron
      hp = 6 * volume / area;
  } // Polyhedron
  else
  {
    log.LogAllError() << "Unsupported cell type in call to HPerpendicular";
    Exit(EXIT_FAILURE);
  }

  return hp;
}

uint64_t
Solver::MapCellLocalNodeIDFromGlobalID(const Cell& cell, uint64_t node_global_id)
{
  size_t imap = 0;
  bool map_found = false;
  for (size_t ai = 0; ai < cell.vertex_ids_.size(); ai++)
  {
    if (node_global_id == cell.vertex_ids_[ai])
    {
      imap = ai;
      map_found = true;
      break;
    }
  }

  if (not map_found) throw std::logic_error(std::string(__FUNCTION__) + ": Mapping failure.");

  return imap;
}

unsigned int
Solver::MapCellFace(const Cell& cur_cell, const Cell& adj_cell, unsigned int f)
{
  const auto& ccface = cur_cell.faces_[f]; // current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids_)
    ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af = 0; af < adj_cell.faces_.size(); af++)
  {
    const auto& acface = adj_cell.faces_[af]; // adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids_)
      acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  } // for adj faces

  if (not map_found) throw std::logic_error(std::string(__FUNCTION__) + ": Mapping failure.");

  return (unsigned int)fmap;
}

} // namespace diffusion
} // namespace opensn
