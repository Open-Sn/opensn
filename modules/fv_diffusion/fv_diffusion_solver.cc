#include "modules/fv_diffusion/fv_diffusion_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/fv_diffusion/fv_diffusion_bndry.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/math/functions/scalar_spatial_material_function.h"

namespace opensn
{
namespace fv_diffusion
{

Solver::Solver(const std::string& in_solver_name)
  : opensn::Solver(in_solver_name, {{"max_iters", int64_t(500)}, {"residual_tolerance", 1.0e-2}})
{
}

fv_diffusion::Solver::~Solver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
}

void
Solver::SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  d_coef_function_ = function;
}

void
Solver::SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  q_ext_function_ = function;
}

void
Solver::SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  sigma_a_function_ = function;
}

// Initialize
void
fv_diffusion::Solver::Initialize()
{
  const std::string fname = "fv_diffusion::Solver::Initialize";
  log.Log() << "\n"
            << program_timer.GetTimeString() << " " << TextName()
            << ": Initializing CFEM Diffusion solver ";

  // Get grid
  grid_ptr_ = GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr_;
  if (grid_ptr_ == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " No grid defined.");

  log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // BIDs
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
          log.Log() << "Boundary " << bndry_name << " set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1, 0.0);
          boundaries_.insert(
            std::make_pair(bndry_id, Boundary{BoundaryType::Dirichlet, {bndry_vals[0], 0., 0.}}));
          log.Log() << "Boundary " << bndry_name << " set to dirichlet.";
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
          log.Log() << "Boundary " << bndry_name << " set to robin." << bndry_vals[0] << ","
                    << bndry_vals[1] << "," << bndry_vals[2];
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries_.insert(
            std::make_pair(bndry_id, Boundary{BoundaryType::Robin, {0.25, 0.5, 0.}}));
          log.Log() << "Boundary " << bndry_name << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Neumann needs 3 values in bndry vals.");
          boundaries_.insert(std::make_pair(
            bndry_id, Boundary{BoundaryType::Robin, {0., bndry_vals[0], bndry_vals[1]}}));
          log.Log() << "Boundary " << bndry_name << " set to neumann." << bndry_vals[0];
          break;
        }
      } // switch boundary type
    }
    else
    {
      boundaries_.insert(std::make_pair(bndry_id, Boundary{BoundaryType::Dirichlet, {0., 0., 0.}}));
      log.Log0Verbose1() << "No boundary preference found for boundary index " << bndry_name
                         << "Dirichlet boundary added with zero boundary value.";
    }
  } // for bndry

  // Make SDM
  sdm_ptr_ = FiniteVolume::New(*grid_ptr_);
  const auto& sdm = *sdm_ptr_;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs_ = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs_ = sdm.GetNumGlobalDOFs(OneDofPerNode);

  log.Log() << "Num local DOFs: " << num_local_dofs_;
  log.Log() << "Num globl DOFs: " << num_globl_dofs_;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs_);
  const auto N = static_cast<int64_t>(num_globl_dofs_);

  A_ = CreateSquareMatrix(n, N);
  x_ = CreateVector(n, N);
  b_ = CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  InitMatrixSparsity(A_, nodal_nnz_in_diag, nodal_nnz_off_diag);

  if (field_functions_.empty())
  {
    std::string solver_name;
    if (not TextName().empty()) solver_name = TextName() + "-";

    std::string text_name = solver_name + "phi";

    auto initial_field_function =
      std::make_shared<FieldFunctionGridBased>(text_name, sdm_ptr_, Unknown(UnknownType::SCALAR));

    field_functions_.push_back(initial_field_function);
    field_function_stack.push_back(initial_field_function);
  } // if not ff set

} // end initialize

// Execute
void
fv_diffusion::Solver::Execute()
{
  log.Log() << "\nExecuting CFEM Diffusion solver";

  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  // Assemble the system
  // P ~ Present cell
  // N ~ Neighbor cell
  log.Log() << "Assembling system: ";
  for (const auto& cell_P : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell_P);
    const double volume_P = cell_mapping.CellVolume(); // Volume of present cell
    const auto& x_cc_P = cell_P.centroid_;

    const auto imat = cell_P.material_id_;

    const double sigma_a = sigma_a_function_->Evaluate(imat, x_cc_P);
    const double q_ext = q_ext_function_->Evaluate(imat, x_cc_P);
    const double D_P = d_coef_function_->Evaluate(imat, x_cc_P);

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

        const double D_N = d_coef_function_->Evaluate(jmat, x_cc_N);

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
  }       // for cell

  log.Log() << "Global assembly";

  MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

  log.Log() << "Done global assembly";

  // Create Krylov Solver
  log.Log() << "Solving: ";
  auto petsc_solver =
    CreateCommonKrylovSolverSetup(A_,
                                  TextName(),
                                  KSPCG,
                                  PCGAMG,
                                  basic_options_("residual_tolerance").FloatValue(),
                                  basic_options_("max_iters").IntegerValue());

  // Solve
  KSPSolve(petsc_solver.ksp, b_, x_);

  UpdateFieldFunctions();

  log.Log() << "Done solving";
}

void
Solver::UpdateFieldFunctions()
{
  auto& ff = *field_functions_.front();
  ff.UpdateFieldVector(x_);
}

} // namespace fv_diffusion
} // namespace opensn
