#include "opensn/modules/FVDiffusion/fv_diffusion_solver.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "opensn/framework/utils/chi_timer.h"

#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "opensn/modules/FVDiffusion/fv_diffusion_bndry.h"

#include "opensn/framework/physics/FieldFunction/fieldfunction_gridbased.h"

#include "opensn/framework/math/SpatialDiscretization/FiniteVolume/FiniteVolume.h"

//============================================= constructor
fv_diffusion::Solver::Solver(const std::string& in_solver_name)
  : chi_physics::Solver(in_solver_name,
                        {{"max_iters", int64_t(500)}, {"residual_tolerance", 1.0e-2}})
{
}

//============================================= destructor
fv_diffusion::Solver::~Solver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
}

//============================================= Initialize
void
fv_diffusion::Solver::Initialize()
{
  const std::string fname = "fv_diffusion::Solver::Initialize";
  Chi::log.Log() << "\n"
                 << Chi::program_timer.GetTimeString() << " " << TextName()
                 << ": Initializing CFEM Diffusion solver ";

  //============================================= Get grid
  grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr_;
  if (grid_ptr_ == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " No grid defined.");

  Chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  //============================================= BIDs
  auto globl_unique_bndry_ids = grid.GetDomainUniqueBoundaryIDs();

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
          boundaries_.insert(
            std::make_pair(bndry_id, Boundary{BoundaryType::Reflecting, {0., 0., 0.}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1, 0.0);
          boundaries_.insert(
            std::make_pair(bndry_id, Boundary{BoundaryType::Dirichlet, {bndry_vals[0], 0., 0.}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to dirichlet.";
          break;
        }
        case BoundaryType::Robin:
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Robin needs 3 values in bndry vals.");
          boundaries_.insert(std::make_pair(
            bndry_id,
            Boundary{BoundaryType::Robin, {bndry_vals[0], bndry_vals[1], bndry_vals[2]}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to robin." << bndry_vals[0] << ","
                         << bndry_vals[1] << "," << bndry_vals[2];
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries_.insert(
            std::make_pair(bndry_id, Boundary{BoundaryType::Robin, {0.25, 0.5, 0.}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Neumann needs 3 values in bndry vals.");
          boundaries_.insert(std::make_pair(
            bndry_id, Boundary{BoundaryType::Robin, {0., bndry_vals[0], bndry_vals[1]}}));
          Chi::log.Log() << "Boundary " << bndry_name << " set to neumann." << bndry_vals[0];
          break;
        }
      } // switch boundary type
    }
    else
    {
      boundaries_.insert(std::make_pair(bndry_id, Boundary{BoundaryType::Dirichlet, {0., 0., 0.}}));
      Chi::log.Log0Verbose1() << "No boundary preference found for boundary index " << bndry_name
                              << "Dirichlet boundary added with zero boundary value.";
    }
  } // for bndry

  //============================================= Make SDM
  sdm_ptr_ = chi_math::spatial_discretization::FiniteVolume::New(*grid_ptr_);
  const auto& sdm = *sdm_ptr_;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs_ = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs_ = sdm.GetNumGlobalDOFs(OneDofPerNode);

  Chi::log.Log() << "Num local DOFs: " << num_local_dofs_;
  Chi::log.Log() << "Num globl DOFs: " << num_globl_dofs_;

  //============================================= Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs_);
  const auto N = static_cast<int64_t>(num_globl_dofs_);

  A_ = chi_math::PETScUtils::CreateSquareMatrix(n, N);
  x_ = chi_math::PETScUtils::CreateVector(n, N);
  b_ = chi_math::PETScUtils::CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  chi_math::PETScUtils::InitMatrixSparsity(A_, nodal_nnz_in_diag, nodal_nnz_off_diag);

  if (field_functions_.empty())
  {
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    using namespace chi_math;
    auto initial_field_function = std::make_shared<chi_physics::FieldFunctionGridBased>(
      text_name,                     // Text name
      sdm_ptr_,                      // Spatial Discretization
      Unknown(UnknownType::SCALAR)); // Unknown/Variable

    field_functions_.push_back(initial_field_function);
    Chi::field_function_stack.push_back(initial_field_function);
  } // if not ff set

} // end initialize

//========================================================== Execute
void
fv_diffusion::Solver::Execute()
{
  Chi::log.Log() << "\nExecuting CFEM Diffusion solver";

  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

#ifdef OPENSN_WITH_LUA
  lua_State* L = Chi::console.GetConsoleState();
#endif

  //============================================= Assemble the system
  // P ~ Present cell
  // N ~ Neighbor cell
  Chi::log.Log() << "Assembling system: ";
  for (const auto& cell_P : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell_P);
    const double volume_P = cell_mapping.CellVolume(); // Volume of present cell
    const auto& x_cc_P = cell_P.centroid_;

    const auto imat = cell_P.material_id_;

#ifdef OPENSN_WITH_LUA
    const double sigma_a = CallLua_iXYZFunction(L, "Sigma_a", imat, x_cc_P);
    const double q_ext = CallLua_iXYZFunction(L, "Q_ext", imat, x_cc_P);
    const double D_P = CallLua_iXYZFunction(L, "D_coef", imat, x_cc_P);

    const int64_t imap = sdm.MapDOF(cell_P, 0);
    MatSetValue(A_, imap, imap, sigma_a * volume_P, ADD_VALUES);
    VecSetValue(b_, imap, q_ext * volume_P, ADD_VALUES);

    for (size_t f = 0; f < cell_P.faces_.size(); ++f)
    {
      const auto& face = cell_P.faces_[f];
      const auto& x_fc = face.centroid_;
      const auto x_PF = x_fc - x_cc_P;
      const auto A_f = cell_mapping.FaceArea(f);
      const auto A_f_n = A_f * face.normal_;

      if (face.has_neighbor_)
      {
        const auto& cell_N = grid.cells[face.neighbor_id_];
        const int jmat = cell_N.material_id_;
        const auto& x_cc_N = cell_N.centroid_;
        const auto x_PN = x_cc_N - x_cc_P;

        const double D_N = CallLua_iXYZFunction(L, "D_coef", jmat, x_cc_N);

        const double w = x_PF.Norm() / x_PN.Norm();
        const double D_f = 1.0 / (w / D_P + (1.0 - w) / D_N);

        const double entry_ii = A_f_n.Dot(D_f * x_PN / x_PN.NormSquare());
        const double entry_ij = -entry_ii;

        const int64_t jmap = sdm.MapDOF(cell_N, 0);
        MatSetValue(A_, imap, imap, entry_ii, ADD_VALUES);
        MatSetValue(A_, imap, jmap, entry_ij, ADD_VALUES);
      } // internal face
      else
      {
        const auto& bndry = boundaries_[face.neighbor_id_];

        if (bndry.type_ == BoundaryType::Robin)
        {
          const auto& aval = bndry.values_[0];
          const auto& bval = bndry.values_[1];
          const auto& fval = bndry.values_[2];

          if (std::fabs(bval) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");

          if (std::fabs(aval) > 1.0e-8) MatSetValue(A_, imap, imap, A_f * aval / bval, ADD_VALUES);
          if (std::fabs(fval) > 1.0e-8) VecSetValue(b_, imap, A_f * fval / bval, ADD_VALUES);
        } // if Robin

        if (bndry.type_ == BoundaryType::Dirichlet)
        {
          const auto& boundary_value = bndry.values_[0];

          const auto& x_cc_N = x_cc_P + 2.0 * x_PF;
          const auto x_PN = x_cc_N - x_cc_P;

          const double D_f = D_P;
          const double entry_ii = A_f_n.Dot(D_f * x_PN / x_PN.NormSquare());

          MatSetValue(A_, imap, imap, entry_ii, ADD_VALUES);
          VecSetValue(b_, imap, entry_ii * boundary_value, ADD_VALUES);
        } // if Dirichlet
      }   // bndry face
    }     // for f
#endif
  }       // for cell

  Chi::log.Log() << "Global assembly";

  MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

  Chi::log.Log() << "Done global assembly";

  //============================================= Create Krylov Solver
  Chi::log.Log() << "Solving: ";
  auto petsc_solver = chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
    A_,                                                // Matrix
    TextName(),                                        // Solver name
    KSPCG,                                             // Solver type
    PCGAMG,                                            // Preconditioner type
    basic_options_("residual_tolerance").FloatValue(), // Relative residual tolerance
    basic_options_("max_iters").IntegerValue()         // Max iterations
  );

  //============================================= Solve
  KSPSolve(petsc_solver.ksp, b_, x_);

  UpdateFieldFunctions();

  Chi::log.Log() << "Done solving";
}
