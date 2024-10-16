// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/diffusion/cfem_diffusion_solver.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/math/functions/scalar_spatial_material_function.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"

namespace opensn
{

OpenSnRegisterObjectAliasInNamespace(diffusion, CFEMSolver, CFEMDiffusionSolver);
OpenSnRegisterSyntaxBlockInNamespace(diffusion,
                                     CFEMBoundaryOptionsBlock,
                                     CFEMDiffusionSolver::BoundaryOptionsBlock);

CFEMDiffusionSolver::CFEMDiffusionSolver(const std::string& name) : DiffusionSolverBase(name)
{
}

InputParameters
CFEMDiffusionSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.AddOptionalParameter<double>("residual_tolerance", 1.0e-2, "Solver relative tolerance");
  params.AddOptionalParameter<int>("max_iters", 500, "Solver relative tolerance");
  return params;
}

InputParameters
CFEMDiffusionSolver::OptionsBlock()
{
  InputParameters params;
  params.AddOptionalParameterArray(
    "boundary_conditions", {}, "An array contain tables for each boundary specification.");
  params.LinkParameterToBlock("boundary_conditions", "CFEMSolver::BoundaryOptionsBlock");
  return params;
}

InputParameters
CFEMDiffusionSolver::BoundaryOptionsBlock()
{
  InputParameters params = DiffusionSolverBase::BoundaryOptionsBlock();
  return params;
}

CFEMDiffusionSolver::CFEMDiffusionSolver(const InputParameters& params)
  : DiffusionSolverBase(params)
{
}

CFEMDiffusionSolver::~CFEMDiffusionSolver()
{
}

void
CFEMDiffusionSolver::SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  d_coef_function_ = function;
}

void
CFEMDiffusionSolver::SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  q_ext_function_ = function;
}

void
CFEMDiffusionSolver::SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  sigma_a_function_ = function;
}

void
CFEMDiffusionSolver::SetOptions(const InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();

  for (size_t p = 0; p < user_params.NumParameters(); ++p)
  {
    const auto& spec = user_params.GetParam(p);
    if (spec.Name() == "boundary_conditions")
    {
      spec.RequireBlockTypeIs(ParameterBlockType::ARRAY);
      for (size_t b = 0; b < spec.NumParameters(); ++b)
      {
        auto bndry_params = BoundaryOptionsBlock();
        bndry_params.AssignParameters(spec.GetParam(b));
        SetBoundaryOptions(bndry_params);
      }
    }
  }
}

void
CFEMDiffusionSolver::SetBoundaryOptions(const InputParameters& params)
{
  const std::string fname = "CFEMSolver::SetBoundaryOptions";

  const auto& user_params = params.ParametersAtAssignment();
  const auto boundary = user_params.GetParamValue<std::string>("boundary");
  const auto bc_type = user_params.GetParamValue<std::string>("type");
  const auto bc_type_lc = LowerCase(bc_type);

  if (bc_type_lc == "reflecting")
  {
    BoundaryInfo bndry_info;
    bndry_info.first = BoundaryType::Reflecting;
    boundary_preferences_.insert(std::make_pair(boundary, bndry_info));
    opensn::log.Log() << "Boundary " << boundary << " set as Reflecting.";
  }
  else if (bc_type_lc == "dirichlet")
  {
    const auto coeffs = user_params.GetParamVectorValue<double>("coeffs");
    if (coeffs.size() < 1)
      throw std::invalid_argument("Expecting one value in the 'coeffs' parameter.");
    auto boundary_value = coeffs[0];
    BoundaryInfo bndry_info;
    bndry_info.first = BoundaryType::Dirichlet;
    bndry_info.second = {boundary_value};
    boundary_preferences_.insert(std::make_pair(boundary, bndry_info));
    opensn::log.Log() << "Boundary " << boundary << " set as Dirichlet with value "
                      << boundary_value;
  }
  else if (bc_type_lc == "neumann")
  {
    const auto coeffs = user_params.GetParamVectorValue<double>("coeffs");
    if (coeffs.size() < 1)
      throw std::invalid_argument("Expecting one value in the 'coeffs' parameter.");
    auto f_value = coeffs[0];
    BoundaryInfo bndry_info;
    bndry_info.first = BoundaryType::Robin;
    bndry_info.second = {0.0, 1.0, f_value};
    boundary_preferences_.insert(std::make_pair(boundary, bndry_info));
    opensn::log.Log() << "Boundary " << boundary << " set as Neumann with D grad(u) dot n = ("
                      << f_value << ") ";
  }
  else if (bc_type_lc == "vacuum")
  {
    BoundaryInfo bndry_info;
    bndry_info.first = BoundaryType::Robin;
    bndry_info.second = {0.25, 0.5, 0.0};
    boundary_preferences_.insert(std::make_pair(boundary, bndry_info));
    opensn::log.Log() << "Boundary " << boundary << " set as Vacuum.";
  }
  else if (bc_type_lc == "robin")
  {
    const auto coeffs = user_params.GetParamVectorValue<double>("coeffs");
    if (coeffs.size() < 3)
      throw std::invalid_argument("Expecting three values in the 'coeffs' parameter.");
    auto a_value = coeffs[0];
    auto b_value = coeffs[1];
    auto f_value = coeffs[2];
    BoundaryInfo bndry_info;
    bndry_info.first = BoundaryType::Robin;
    bndry_info.second = {a_value, b_value, f_value};
    boundary_preferences_.insert(std::make_pair(boundary, bndry_info));
    opensn::log.Log() << "Boundary " << boundary << " set as Robin with a,b,f = (" << a_value << ","
                      << b_value << "," << f_value << ") ";
  }
  else
    throw std::invalid_argument(fname + ": Unsupported boundary type '" + bc_type + "'.");
}

void
CFEMDiffusionSolver::Initialize()
{
  const std::string fname = "CFEMSolver::Initialize";
  log.Log() << "\n"
            << program_timer.GetTimeString() << " " << TextName()
            << ": Initializing CFEM Diffusion solver ";

  // Get grid
  grid_ptr_ = GetCurrentMesh();
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
            std::make_pair(bndry_id, Boundary{BoundaryType::Reflecting, {0.0, 0.0, 0.0}}));
          log.Log() << "Boundary " << bndry_name << " set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty())
            bndry_vals.resize(1, 0.0);
          boundaries_.insert(
            std::make_pair(bndry_id, Boundary{BoundaryType::Dirichlet, {bndry_vals[0], 0.0, 0.0}}));
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
            bndry_id, Boundary{BoundaryType::Robin, {0.0, bndry_vals[0], bndry_vals[1]}}));
          log.Log() << "Boundary " << bndry_name << " set to neumann." << bndry_vals[0];
          break;
        }
      } // switch boundary type
    }
    else
    {
      boundaries_.insert(
        std::make_pair(bndry_id, Boundary{BoundaryType::Dirichlet, {0.0, 0.0, 0.0}}));
      log.Log0Verbose1() << "No boundary preference found for boundary index " << bndry_id
                         << "Dirichlet boundary added with zero boundary value.";
    }
  } // for bndry

  // Make SDM
  sdm_ptr_ = PieceWiseLinearContinuous::New(*grid_ptr_);
  const auto& sdm = *sdm_ptr_;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs_ = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_global_dofs_ = sdm.GetNumGlobalDOFs(OneDofPerNode);

  log.Log() << "Num local DOFs: " << num_local_dofs_;
  log.Log() << "Num globl DOFs: " << num_global_dofs_;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs_);
  const auto N = static_cast<int64_t>(num_global_dofs_);

  A_ = CreateSquareMatrix(n, N);
  x_ = CreateVector(n, N);
  b_ = CreateVector(n, N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  InitMatrixSparsity(A_, nodal_nnz_in_diag, nodal_nnz_off_diag);

  InitFieldFunctions();
}

void
CFEMDiffusionSolver::Execute()
{
  log.Log() << "\nExecuting CFEM Diffusion solver";

  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  // Assemble the system
  log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    const auto imat = cell.material_id;
    const size_t num_nodes = cell_mapping.NumNodes();
    MatDbl Acell(num_nodes, std::vector<double>(num_nodes, 0.0));
    std::vector<double> cell_rhs(num_nodes, 0.0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
      {
        double entry_aij = 0.0;
        for (size_t qp : fe_vol_data.QuadraturePointIndices())
        {
          entry_aij += (d_coef_function_->Evaluate(imat, fe_vol_data.QPointXYZ(qp)) *
                          fe_vol_data.ShapeGrad(i, qp).Dot(fe_vol_data.ShapeGrad(j, qp)) +
                        sigma_a_function_->Evaluate(imat, fe_vol_data.QPointXYZ(qp)) *
                          fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp)) *
                       fe_vol_data.JxW(qp);
        } // for qp
        Acell[i][j] = entry_aij;
      } // for j
      for (size_t qp : fe_vol_data.QuadraturePointIndices())
        cell_rhs[i] += q_ext_function_->Evaluate(imat, fe_vol_data.QPointXYZ(qp)) *
                       fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
    } // for i

    // Flag nodes for being on a boundary
    std::vector<int> dirichlet_count(num_nodes, 0);
    std::vector<double> dirichlet_value(num_nodes, 0.0);

    const size_t num_faces = cell.faces.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      // not a boundary face
      if (face.has_neighbor)
        continue;

      const auto& bndry = boundaries_[face.neighbor_id];

      // Robin boundary
      if (bndry.type == BoundaryType::Robin)
      {
        const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
        const size_t num_face_nodes = face.vertex_ids.size();

        const auto& aval = bndry.values[0];
        const auto& bval = bndry.values[1];
        const auto& fval = bndry.values[2];

        log.Log0Verbose1() << "Boundary  set as Robin with a,b,f = (" << aval << "," << bval << ","
                           << fval << ") ";
        // true Robin when a!=0, otherwise, it is a Neumann:
        // Assert if b=0
        if (std::fabs(bval) < 1e-8)
          throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");

        // loop over nodes of that face
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const uint i = cell_mapping.MapFaceNode(f, fi);

          double entry_rhsi = 0.0;
          for (size_t qp : fe_srf_data.QuadraturePointIndices())
            entry_rhsi += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.JxW(qp);
          cell_rhs[i] += fval / bval * entry_rhsi;

          // only do this part if true Robin (i.e., a!=0)
          if (std::fabs(aval) > 1.0e-8)
          {
            for (size_t fj = 0; fj < num_face_nodes; ++fj)
            {
              const uint j = cell_mapping.MapFaceNode(f, fj);

              double entry_aij = 0.0;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                entry_aij += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(j, qp) *
                             fe_srf_data.JxW(qp);
              Acell[i][j] += aval / bval * entry_aij;
            } // for fj
          }   // end true Robin
        }     // for fi
      }       // if Robin

      // Dirichlet boundary
      if (bndry.type == BoundaryType::Dirichlet)
      {
        const size_t num_face_nodes = face.vertex_ids.size();

        const auto& boundary_value = bndry.values[0];

        // loop over nodes of that face
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const uint i = cell_mapping.MapFaceNode(f, fi);
          dirichlet_count[i] += 1;
          dirichlet_value[i] += boundary_value;
        } // for fi
      }   // if Dirichlet

    } // for face f

    // Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    // Assembly into system
    for (size_t i = 0; i < num_nodes; ++i)
    {
      if (dirichlet_count[i] > 0) // if Dirichlet boundary node
      {
        MatSetValue(A_, imap[i], imap[i], 1.0, ADD_VALUES);
        // because we use CFEM, a given node is common to several faces
        const double aux = dirichlet_value[i] / dirichlet_count[i];
        VecSetValue(b_, imap[i], aux, ADD_VALUES);
      }
      else
      {
        for (size_t j = 0; j < num_nodes; ++j)
        {
          if (dirichlet_count[j] == 0) // not related to a dirichlet node
            MatSetValue(A_, imap[i], imap[j], Acell[i][j], ADD_VALUES);
          else
          {
            const double aux = dirichlet_value[j] / dirichlet_count[j];
            cell_rhs[i] -= Acell[i][j] * aux;
          }
        } // for j
        VecSetValue(b_, imap[i], cell_rhs[i], ADD_VALUES);
      }
    } // for i
  }   // for cell

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
                                  0.0,
                                  basic_options_("residual_tolerance").FloatValue(),
                                  basic_options_("max_iters").IntegerValue());

  // Solve
  KSPSolve(petsc_solver.ksp, b_, x_);

  UpdateFieldFunctions();

  log.Log() << "Done solving";
}

} // namespace opensn
