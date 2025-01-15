// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/diffusion/dfem_diffusion_solver.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/functions/scalar_spatial_material_function.h"

namespace opensn
{

OpenSnRegisterObjectAliasInNamespace(diffusion, DFEMSolver, DFEMDiffusionSolver);

std::shared_ptr<DFEMDiffusionSolver>
DFEMDiffusionSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<DFEMDiffusionSolver>("diffusion::DFEMDiffusionSolver", params);
}

DFEMDiffusionSolver::DFEMDiffusionSolver(const std::string& name) : DiffusionSolverBase(name)
{
}

InputParameters
DFEMDiffusionSolver::GetInputParameters()
{
  InputParameters params = DiffusionSolverBase::GetInputParameters();
  return params;
}

InputParameters
DFEMDiffusionSolver::OptionsBlock()
{
  InputParameters params = DiffusionSolverBase::OptionsBlock();
  return params;
}

InputParameters
DFEMDiffusionSolver::BoundaryOptionsBlock()
{
  InputParameters params = DiffusionSolverBase::BoundaryOptionsBlock();
  return params;
}

DFEMDiffusionSolver::DFEMDiffusionSolver(const InputParameters& params)
  : DiffusionSolverBase(params)
{
}

DFEMDiffusionSolver::~DFEMDiffusionSolver()
{
}

void
DFEMDiffusionSolver::SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  d_coef_function_ = function;
}

void
DFEMDiffusionSolver::SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  q_ext_function_ = function;
}

void
DFEMDiffusionSolver::SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  sigma_a_function_ = function;
}

void
DFEMDiffusionSolver::SetOptions(const InputParameters& params)
{
  for (size_t p = 0; p < params.NumParameters(); ++p)
  {
    const auto& spec = params.GetParam(p);
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
DFEMDiffusionSolver::SetBoundaryOptions(const InputParameters& params)
{
  const std::string fname = "DFEMSolver::SetBoundaryOptions";

  const auto boundary = params.GetParamValue<std::string>("boundary");
  const auto bc_type = params.GetParamValue<std::string>("type");
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
    const auto coeffs = params.GetParamVectorValue<double>("coeffs");
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
    const auto coeffs = params.GetParamVectorValue<double>("coeffs");
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
    const auto coeffs = params.GetParamVectorValue<double>("coeffs");
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
DFEMDiffusionSolver::Initialize()
{
  const std::string fname = "DFEMSolver::Initialize";
  log.Log() << "\n"
            << program_timer.GetTimeString() << " " << Name()
            << ": Initializing DFEM Diffusion solver ";

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
      log.Log0Verbose1() << "No boundary preference found for boundary index " << bndry_name
                         << "Dirichlet boundary added with zero boundary value.";
    }
  } // for bndry

  // Make SDM
  sdm_ptr_ = PieceWiseLinearDiscontinuous::New(*grid_ptr_);
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
DFEMDiffusionSolver::Execute()
{
  log.Log() << "\nExecuting DFEM IP Diffusion solver";

  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  // Assemble the system
  // is this needed?
  VecSet(b_, 0.0);

  log.Log() << "Assembling system: ";

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto cc_nodes = cell_mapping.GetNodeLocations();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    const auto imat = cell.material_id;

    // Assemble volumetric terms
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOF(cell, i);

      for (size_t j = 0; j < num_nodes; ++j)
      {
        const int64_t jmap = sdm.MapDOF(cell, j);
        double entry_aij = 0.0;
        for (size_t qp : fe_vol_data.QuadraturePointIndices())
        {
          entry_aij += (d_coef_function_->Evaluate(imat, fe_vol_data.QPointXYZ(qp)) *
                          fe_vol_data.ShapeGrad(i, qp).Dot(fe_vol_data.ShapeGrad(j, qp)) +
                        sigma_a_function_->Evaluate(imat, fe_vol_data.QPointXYZ(qp)) *
                          fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp)) *
                       fe_vol_data.JxW(qp);
        } // for qp
        MatSetValue(A_, imap, jmap, entry_aij, ADD_VALUES);
      } // for j
      double entry_rhs_i = 0.0;
      for (size_t qp : fe_vol_data.QuadraturePointIndices())
        entry_rhs_i += q_ext_function_->Evaluate(imat, fe_vol_data.QPointXYZ(qp)) *
                       fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
      VecSetValue(b_, imap, entry_rhs_i, ADD_VALUES);
    } // for i

    // Assemble face terms
    const size_t num_faces = cell.faces.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& n_f = face.normal;
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);

      const double hm = HPerpendicular(cell, f);

      // interior face
      if (face.has_neighbor)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id];
        const auto& adj_cell_mapping = sdm.GetCellMapping(adj_cell);
        const auto ac_nodes = adj_cell_mapping.GetNodeLocations();
        const size_t acf = MeshContinuum::MapCellFace(cell, adj_cell, f);
        const double hp_neigh = HPerpendicular(adj_cell, acf);

        const auto imat_neigh = adj_cell.material_id;

        // Compute Ckappa IP
        double Ckappa = 1.0;
        if (cell.Type() == CellType::SLAB)
          Ckappa = 2.0;
        if (cell.Type() == CellType::POLYGON)
          Ckappa = 2.0;
        if (cell.Type() == CellType::POLYHEDRON)
          Ckappa = 4.0;

        // Assembly penalty terms
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const int64_t imap = sdm.MapDOF(cell, i);

          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
            const int64_t jmmap = sdm.MapDOF(cell, jm);

            const int jp =
              MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fj); // j-plus
            const int64_t jpmap = sdm.MapDOF(adj_cell, jp);

            double aij = 0.0;
            for (size_t qp : fe_srf_data.QuadraturePointIndices())
              aij +=
                Ckappa *
                (d_coef_function_->Evaluate(imat, fe_srf_data.QPointXYZ(qp)) / hm +
                 d_coef_function_->Evaluate(imat_neigh, fe_srf_data.QPointXYZ(qp)) / hp_neigh) /
                2.0 * fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                fe_srf_data.JxW(qp);

            MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
          } // for fj
        }   // for fi

        // Assemble gradient terms
        // For the following comments we use the notation:
        // Dk = 0.5* n dot nabla bk

        // {{D d_n b_i}}[[Phi]]
        // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-

        // loop over node of current cell (gradient of b_i)
        for (int i = 0; i < num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOF(cell, i);

          // loop over faces
          for (int fj = 0; fj < num_face_nodes; ++fj)
          {
            const int jm = cell_mapping.MapFaceNode(f, fj); // j-minus
            const int64_t jmmap = sdm.MapDOF(cell, jm);
            const int jp =
              MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fj); // j-plus
            const int64_t jpmap = sdm.MapDOF(adj_cell, jp);

            Vector3 vec_aij;
            for (size_t qp : fe_srf_data.QuadraturePointIndices())
              vec_aij += d_coef_function_->Evaluate(imat, fe_srf_data.QPointXYZ(qp)) *
                         fe_srf_data.ShapeValue(jm, qp) * fe_srf_data.ShapeGrad(i, qp) *
                         fe_srf_data.JxW(qp);
            const double aij = -0.5 * n_f.Dot(vec_aij);

            MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
            MatSetValue(A_, imap, jpmap, -aij, ADD_VALUES);
          } // for fj
        }   // for i

        // {{D d_n Phi}}[[b_i]]
        // 0.5*D* n dot (b_i^+ - b_i^-)*nabla b_j^-
        for (int fi = 0; fi < num_face_nodes; ++fi)
        {
          const int im = cell_mapping.MapFaceNode(f, fi); // i-minus
          const int64_t immap = sdm.MapDOF(cell, im);

          const int ip = MapFaceNodeDisc(cell, adj_cell, cc_nodes, ac_nodes, f, acf, fi); // i-plus
          const int64_t ipmap = sdm.MapDOF(adj_cell, ip);

          for (int j = 0; j < num_nodes; ++j)
          {
            const int64_t jmap = sdm.MapDOF(cell, j);

            Vector3 vec_aij;
            for (size_t qp : fe_srf_data.QuadraturePointIndices())
              vec_aij += d_coef_function_->Evaluate(imat, fe_srf_data.QPointXYZ(qp)) *
                         fe_srf_data.ShapeValue(im, qp) * fe_srf_data.ShapeGrad(j, qp) *
                         fe_srf_data.JxW(qp);
            const double aij = -0.5 * n_f.Dot(vec_aij);

            MatSetValue(A_, immap, jmap, aij, ADD_VALUES);
            MatSetValue(A_, ipmap, jmap, -aij, ADD_VALUES);
          } // for j
        }   // for fi

      } // internal face
      else
      { // boundary face
        const auto& bndry = boundaries_[face.neighbor_id];
        // Robin boundary
        if (bndry.type == BoundaryType::Robin)
        {
          const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);

          const auto& aval = bndry.values[0];
          const auto& bval = bndry.values[1];
          const auto& fval = bndry.values[2];

          log.Log0Verbose1() << "Boundary  set as Robin with a,b,f = (" << aval << "," << bval
                             << "," << fval << ") ";
          // true Robin when a!=0, otherwise, it is a Neumann:
          // Assert if b=0
          if (std::fabs(bval) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");

          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const uint i = cell_mapping.MapFaceNode(f, fi);
            const int64_t ir = sdm.MapDOF(cell, i);

            if (std::fabs(aval) >= 1.0e-12)
            {
              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const uint j = cell_mapping.MapFaceNode(f, fj);
                const int64_t jr = sdm.MapDOF(cell, j);

                double aij = 0.0;
                for (size_t qp : fe_srf_data.QuadraturePointIndices())
                  aij += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(j, qp) *
                         fe_srf_data.JxW(qp);
                aij *= (aval / bval);

                MatSetValue(A_, ir, jr, aij, ADD_VALUES);
              } // for fj
            }   // if a nonzero

            if (std::fabs(fval) >= 1.0e-12)
            {
              double rhs_val = 0.0;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                rhs_val += fe_srf_data.ShapeValue(i, qp) * fe_srf_data.JxW(qp);
              rhs_val *= (fval / bval);

              VecSetValue(b_, ir, rhs_val, ADD_VALUES);
            } // if f nonzero
          }   // for fi
        }     // Robin BC
        else if (bndry.type == BoundaryType::Dirichlet)
        {
          const double bc_value = bndry.values[0];
          // Compute kappa
          double Ckappa = 2.0;
          if (cell.Type() == CellType::SLAB)
            Ckappa = 4.0; // fmax(4.0*Dg/hm,0.25);
          if (cell.Type() == CellType::POLYGON)
            Ckappa = 4.0;
          if (cell.Type() == CellType::POLYHEDRON)
            Ckappa = 8.0;

          // Assembly penalty terms
          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const uint i = cell_mapping.MapFaceNode(f, fi);
            const int64_t imap = sdm.MapDOF(cell, i);

            for (size_t fj = 0; fj < num_face_nodes; ++fj)
            {
              const uint jm = cell_mapping.MapFaceNode(f, fj);
              const int64_t jmmap = sdm.MapDOF(cell, jm);

              double aij = 0.0;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                aij += Ckappa * d_coef_function_->Evaluate(imat, fe_srf_data.QPointXYZ(qp)) / hm *
                       fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(jm, qp) *
                       fe_srf_data.JxW(qp);
              double aij_bc_value = aij * bc_value;

              MatSetValue(A_, imap, jmmap, aij, ADD_VALUES);
              VecSetValue(b_, imap, aij_bc_value, ADD_VALUES);
            } // for fj
          }   // for fi

          // Assemble gradient terms
          // For the following comments we use the notation:
          // Dk = 0.5* n dot nabla bk

          // 0.5*D* n dot (b_j^+ - b_j^-)*nabla b_i^-
          for (size_t i = 0; i < num_nodes; ++i)
          {
            const int64_t imap = sdm.MapDOF(cell, i);

            for (size_t j = 0; j < num_nodes; ++j)
            {
              const int64_t jmap = sdm.MapDOF(cell, j);

              Vector3 vec_aij;
              for (size_t qp : fe_srf_data.QuadraturePointIndices())
                vec_aij += (fe_srf_data.ShapeValue(j, qp) * fe_srf_data.ShapeGrad(i, qp) +
                            fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeGrad(j, qp)) *
                           fe_srf_data.JxW(qp) *
                           d_coef_function_->Evaluate(imat, fe_srf_data.QPointXYZ(qp));

              const double aij = -n_f.Dot(vec_aij);
              double aij_bc_value = aij * bc_value;

              MatSetValue(A_, imap, jmap, aij, ADD_VALUES);
              VecSetValue(b_, imap, aij_bc_value, ADD_VALUES);
            }   // for fj
          }     // for i
        }       // Dirichlet BC
        else {} // else BC
      }         // boundary face
    }           // for face f
  }             // for cell

  log.Log() << "Global assembly";

  MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);

  //  MatView(A, PETSC_VIEWER_STDERR_WORLD);
  //
  //  PetscViewer viewer;
  //  PetscViewerASCIIOpen(opensn::mpi_comm,"A.m",&viewer);
  //  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  //  MatView(A,viewer);
  //  PetscViewerPopFormat(viewer);
  //  PetscViewerDestroy(&viewer);

  log.Log() << "Done global assembly";

  // Create Krylov Solver
  log.Log() << "Solving: ";
  auto petsc_solver =
    CreateCommonKrylovSolverSetup(A_,
                                  Name(),
                                  KSPCG,
                                  PCGAMG,
                                  0.0,
                                  basic_options_("residual_tolerance").FloatValue(),
                                  basic_options_("max_iters").IntegerValue());

  // Solve
  KSPSolve(petsc_solver.ksp, b_, x_);

  log.Log() << "Done solving";

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  sdm.LocalizePETScVector(x_, field_, OneDofPerNode);

  field_functions_.front()->UpdateFieldVector(field_);
}

double
DFEMDiffusionSolver::HPerpendicular(const Cell& cell, unsigned int f)
{
  const auto& sdm = *sdm_ptr_;

  const auto& cell_mapping = sdm.GetCellMapping(cell);
  double hp;

  const size_t num_faces = cell.faces.size();
  const size_t num_vertices = cell.vertex_ids.size();

  const double volume = cell_mapping.CellVolume();
  const double face_area = cell_mapping.FaceArea(f);

  /**Lambda to compute surface area.*/
  auto ComputeSurfaceArea = [&cell_mapping, &num_faces]()
  {
    double surface_area = 0.0;
    for (size_t fr = 0; fr < num_faces; ++fr)
      surface_area += cell_mapping.FaceArea(fr);

    return surface_area;
  };

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == CellType::SLAB)
    hp = volume / 2.0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == CellType::POLYGON)
  {
    if (num_faces == 3)
      hp = 2.0 * volume / face_area;
    else if (num_faces == 4)
      hp = volume / face_area;
    else // Nv > 4
    {
      const double surface_area = ComputeSurfaceArea();

      if (num_faces % 2 == 0)
        hp = 4.0 * volume / surface_area;
      else
      {
        hp = 2.0 * volume / surface_area;
        hp +=
          sqrt(2.0 * volume /
               (static_cast<double>(num_faces) * sin(2.0 * M_PI / static_cast<double>(num_faces))));
      }
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == CellType::POLYHEDRON)
  {
    const double surface_area = ComputeSurfaceArea();

    if (num_faces == 4) // Tet
      hp = 3 * volume / surface_area;
    else if (num_faces == 6 and num_vertices == 8) // Hex
      hp = volume / surface_area;
    else // Polyhedron
      hp = 6 * volume / surface_area;
  } // Polyhedron
  else
    throw std::logic_error(
      "DFEMSolver::HPerpendicular: Unsupported cell type in call to HPerpendicular");

  return hp;
}

int
DFEMDiffusionSolver::MapFaceNodeDisc(const Cell& cur_cell,
                                     const Cell& adj_cell,
                                     const std::vector<Vector3>& cc_node_locs,
                                     const std::vector<Vector3>& ac_node_locs,
                                     size_t ccf,
                                     size_t acf,
                                     size_t ccfi,
                                     double epsilon)
{
  const auto& sdm = *sdm_ptr_;

  const auto& cur_cell_mapping = sdm.GetCellMapping(cur_cell);
  const auto& adj_cell_mapping = sdm.GetCellMapping(adj_cell);

  const int i = cur_cell_mapping.MapFaceNode(ccf, ccfi);
  const auto& node_i_loc = cc_node_locs[i];

  const size_t adj_face_num_nodes = adj_cell_mapping.NumFaceNodes(acf);

  for (size_t fj = 0; fj < adj_face_num_nodes; ++fj)
  {
    const int j = adj_cell_mapping.MapFaceNode(acf, fj);
    if ((node_i_loc - ac_node_locs[j]).NormSquare() < epsilon)
      return j;
  }

  throw std::logic_error("DFEMSolver::MapFaceNodeDisc: Mapping failure.");
}

} // namespace opensn
