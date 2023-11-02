#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"

//============================================= constructor
mg_diffusion::Solver::Solver(const std::string& in_solver_name)
  : chi_physics::Solver(in_solver_name,
                        {{"max_inner_iters", int64_t(500)},
                         {"residual_tolerance", 1.0e-2},
                         {"verbose_level", int64_t(0)},
                         {"thermal_flux_tolerance", 1.0e-2},
                         {"max_thermal_iters", int64_t(500)},
                         {"do_two_grid", false}})
{
}

//============================================= destructor
mg_diffusion::Solver::~Solver()
{
  for (uint g = 0; g < num_groups_; ++g)
  {
    VecDestroy(&x_[g]);
    VecDestroy(&bext_[g]);
    MatDestroy(&A_[g]);
  }
  VecDestroy(&b_);

  if (last_fast_group_ < num_groups_)
  {
    VecDestroy(&thermal_dphi_);
    for (uint g = last_fast_group_; g < num_groups_; ++g)
      VecDestroy(&x_old_[g]);
  }

  if (do_two_grid_)
  {
    VecDestroy(&x_[num_groups_]);
    MatDestroy(&A_[num_groups_]);
  }
}

//============================================= Initialize
void
mg_diffusion::Solver::Initialize()
{
  Chi::log.Log() << "\n"
                 << Chi::program_timer.GetTimeString() << " " << TextName()
                 << ": Initializing CFEM Multigroup Diffusion solver ";

  //============================================= Get grid
  grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr_;
  if (grid_ptr_ == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " No grid defined.");

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //================================================== Add unique material ids
  std::set<int> unique_material_ids;
  int invalid_mat_cell_count = 0;
  for (auto& cell : grid.local_cells)
  {
    unique_material_ids.insert(cell.material_id_);
    if (cell.material_id_ < 0) ++invalid_mat_cell_count;
  }
  const auto& ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t cell_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[cell_id];
    unique_material_ids.insert(cell.material_id_);
    if (cell.material_id_ < 0) ++invalid_mat_cell_count;
  }

  if (invalid_mat_cell_count > 0)
  {
    Chi::log.LogAllWarning() << "Number of invalid material cells: " << invalid_mat_cell_count;
  }

  //================================================== Initialize materials
  mg_diffusion::Solver::Initialize_Materials(unique_material_ids);

  //============================================= BIDs
  auto globl_unique_bndry_ids = grid.GetDomainUniqueBoundaryIDs();
  mg_diffusion::Solver::Set_BCs(globl_unique_bndry_ids);

  //============================================= Make SDM
  sdm_ptr_ = chi_math::spatial_discretization::PieceWiseLinearContinuous::New(*grid_ptr_);
  const auto& sdm = *sdm_ptr_;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs_ = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs_ = sdm.GetNumGlobalDOFs(OneDofPerNode);

  Chi::log.Log() << "Num local DOFs: " << num_local_dofs_;
  Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs_;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs_);
  const auto N = static_cast<int64_t>(num_globl_dofs_);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  unsigned int i_two_grid = do_two_grid_ ? 1 : 0;
  //  std::cout << "i_two_grid = " << i_two_grid << std::endl;

  A_.resize(num_groups_ + i_two_grid, nullptr);
  x_.resize(num_groups_ + i_two_grid, nullptr);
  bext_.resize(num_groups_, nullptr);

  auto ghost_dof_indices = sdm.GetGhostDOFIndices(OneDofPerNode);

  for (uint g = 0; g < num_groups_; ++g)
  {
    // x[g] = chi_math::PETScUtils::CreateVector(n,N);
    x_[g] = chi_math::PETScUtils::CreateVectorWithGhosts(
      n, N, static_cast<int64_t>(ghost_dof_indices.size()), ghost_dof_indices);
    VecSet(x_[g], 0.0);
    bext_[g] = chi_math::PETScUtils::CreateVector(n, N);

    A_[g] = chi_math::PETScUtils::CreateSquareMatrix(n, N);
    chi_math::PETScUtils::InitMatrixSparsity(A_[g], nodal_nnz_in_diag, nodal_nnz_off_diag);
  }
  // initialize b
  VecDuplicate(bext_.front(), &b_);
  // initialize old flux for thermal groups
  if (last_fast_group_ < num_groups_)
  {
    VecDuplicate(x_.front(), &thermal_dphi_);
    x_old_.resize(num_groups_, nullptr);
    for (uint g = 0; g < num_groups_; ++g)
    {
      VecDuplicate(x_.front(), &x_old_[g]); // jcr is x_old like bext or like x?
      VecSet(x_old_[g], 0.0);
    }
  }
  // add two-grid mat and vec, if needed
  if (do_two_grid_)
  {
    A_[num_groups_] = chi_math::PETScUtils::CreateSquareMatrix(n, N);
    chi_math::PETScUtils::InitMatrixSparsity(
      A_[num_groups_], nodal_nnz_in_diag, nodal_nnz_off_diag);
    VecDuplicate(x_.front(), &x_[num_groups_]); // jcr needed?
    //    x[num_groups] = chi_math::PETScUtils::CreateVectorWithGhosts(n,N,
    //                                                        static_cast<int64_t>(ghost_dof_indices.size()),
    //                                                        ghost_dof_indices);
  }

  if (do_two_grid_) mg_diffusion::Solver::Compute_TwoGrid_VolumeFractions();

  //============================================= Create Mats and ExtVecs
  mg_diffusion::Solver::Assemble_A_bext();

  //============================================= Volume fraction for two-grid
  //                                              update

  //============================================= Field Function
  if (field_functions_.empty())
  {
    for (uint g = 0; g < mg_diffusion::Solver::num_groups_; ++g)
    {
      std::string solver_name;
      if (not TextName().empty()) solver_name = TextName() + "-";

      char buff[100];
      int dummy = snprintf(buff, 4, "%03d", g);

      std::string text_name = solver_name + "Flux_g" + std::string(buff);

      using namespace chi_math;
      auto initial_field_function = std::make_shared<chi_physics::FieldFunctionGridBased>(
        text_name,                     // Text name
        sdm_ptr_,                      // Spatial Discretization
        Unknown(UnknownType::SCALAR)); // Unknown Manager

      field_functions_.push_back(initial_field_function);
      Chi::field_function_stack.push_back(initial_field_function);
    } // for g
  }   // if not ff set

} // end initialize
