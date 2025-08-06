// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#include "framework/runtime.h"
#include "framework/data_types/vector_ghost_communicator/vector_ghost_communicator.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "modules/diffusion/diffusion.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/lbs_keigen_acceleration.h"

namespace opensn
{
InputParameters
LBSKEigenAcceleration::GetInputParameters()
{
  InputParameters params;

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the acceleration solver. This name will be used in status "
    "messages and verbose iterative convergence monitors.");

  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "An existing discrete ordinates problem");
  params.AddOptionalParameter("l_abs_tol", 1.0e-10, "Absolute residual tolerance");
  params.AddOptionalParameter("max_iters", 100, "Maximum allowable iterations");
  params.AddOptionalParameter("verbose", false, "If true, enables verbose output");
  params.AddOptionalParameter("petsc_options", std::string("ssss"), "Additional PETSc options");
  params.AddOptionalParameter(
    "pi_max_its", 50, "Maximum allowable iterations for inner power iterations");
  params.AddOptionalParameter(
    "pi_k_tol", 1.0e-10, "k-eigenvalue tolerance for the inner power iterations");

  params.ConstrainParameterRange("l_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("max_iters", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("pi_max_its", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("pi_k_tol", AllowableRangeLowLimit::New(1.0e-18));

  return params;
}

LBSKEigenAcceleration::LBSKEigenAcceleration(const InputParameters& params)
  : do_problem_(*params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem")),
    l_abs_tol_(params.GetParamValue<double>("l_abs_tol")),
    max_iters_(params.GetParamValue<int>("max_iters")),
    verbose_(params.GetParamValue<bool>("verbose")),
    petsc_options_(params.GetParamValue<std::string>("petsc_options")),
    pi_max_its_(params.GetParamValue<int>("pi_max_its")),
    pi_k_tol_(params.GetParamValue<double>("pi_k_tol")),
    groupsets_(do_problem_.GetGroupsets()),
    front_gs_(groupsets_.front()),
    q_moments_local_(do_problem_.GetQMomentsLocal()),
    phi_old_local_(do_problem_.GetPhiOldLocal()),
    phi_new_local_(do_problem_.GetPhiNewLocal()),
    name_(params.GetParamValue<std::string>("name"))
{
  if (do_problem_.GetGroupsets().size() != 1)
    throw std::logic_error(
      "LBS acceleration schemes are only implemented for problems with a single groupset.");

  // If using the AAH solver with one sweep, a few iterations need to be done
  // to get rid of the junk in the unconverged lagged angular fluxes.  Five
  // sweeps is a guess at how many initial sweeps are necessary.
  if (const auto do_problem = dynamic_cast<DiscreteOrdinatesProblem*>(&do_problem_))
    if (do_problem->GetSweepType() == "AAH" and front_gs_.max_iterations == 1)
      throw std::logic_error("The AAH solver is not stable for single-sweep methods due to "
                             "the presence of lagged angular fluxes.  Multiple sweeps are "
                             "allowed, however, the number of sweeps required to get sensible "
                             "results is not well studied and problem dependent.");
}

void
LBSKEigenAcceleration::Initialize(PowerIterationKEigenSolver& solver)
{
  solver_ = &solver;
  Initialize();
}

std::vector<int64_t>
LBSKEigenAcceleration::MakePWLDGhostIndices(const SpatialDiscretization& pwld,
                                            const UnknownManager& uk_man)
{
  std::set<int64_t> ghost_ids;
  const auto& grid = *pwld.GetGrid();
  for (const uint64_t ghost_id : grid.cells.GetGhostGlobalIDs())
  {
    const auto& cell = grid.cells[ghost_id];
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    for (int i = 0; i < cell_mapping.GetNumNodes(); ++i)
      for (int u = 0; u < uk_man.GetNumberOfUnknowns(); ++u)
        for (int c = 0; c < uk_man.unknowns[u].GetNumComponents(); ++c)
          ghost_ids.insert(pwld.MapDOF(cell, i, uk_man, u, c));
  }
  return {ghost_ids.begin(), ghost_ids.end()};
}

LBSKEigenAcceleration::GhostInfo
LBSKEigenAcceleration::MakePWLDGhostInfo(const SpatialDiscretization& pwld,
                                         const UnknownManager& uk_man)

{
  const auto num_local_dofs = pwld.GetNumLocalDOFs(uk_man);
  const auto num_global_dofs = pwld.GetNumGlobalDOFs(uk_man);

  // Build list of global IDs
  auto ghost_ids = MakePWLDGhostIndices(pwld, uk_man);

  // Create the ghost communicator
  auto vec_ghost_comm = std::make_shared<VectorGhostCommunicator>(
    num_local_dofs, num_global_dofs, ghost_ids, opensn::mpi_comm);

  // Create the mapping
  int64_t k = 0;
  std::map<int64_t, int64_t> ghost_global_to_local_map;
  for (const int64_t ghost_id : ghost_ids)
    ghost_global_to_local_map[ghost_id] = static_cast<int64_t>(num_local_dofs + k++);

  return {vec_ghost_comm, ghost_global_to_local_map};
}

void
LBSKEigenAcceleration::InitializeLinearContinuous()
{
  const auto& sdm = do_problem_.GetSpatialDiscretization();
  pwlc_ptr_ = PieceWiseLinearContinuous::New(sdm.GetGrid());
  ghost_info_ = MakePWLDGhostInfo(*pwlc_ptr_, do_problem_.GetUnknownManager());
}

void
LBSKEigenAcceleration::NodallyAveragedPWLDVector(const std::vector<double>& input,
                                                 std::vector<double>& output) const
{
  const auto& pwld_sdm = do_problem_.GetSpatialDiscretization();
  const auto& pwlc_sdm = diffusion_solver_->GetSpatialDiscretization();
  const auto& uk_man = do_problem_.GetUnknownManager();

  const auto& vgc = ghost_info_.vector_ghost_communicator;
  const auto& dfem_dof_global2local_map = ghost_info_.ghost_global_id_2_local_map;

  auto input_with_ghosts = vgc->MakeGhostedVector(input);
  vgc->CommunicateGhostEntries(input_with_ghosts);

  const auto& grid = pwld_sdm.GetGrid();

  const size_t num_unknowns = uk_man.unknowns.size();

  const size_t num_cfem_local_dofs = pwlc_sdm.GetNumLocalAndGhostDOFs(uk_man);

  std::vector<double> cont_input(num_cfem_local_dofs, 0.0);
  std::vector<double> cont_input_ctr(num_cfem_local_dofs, 0.0);

  std::map<int64_t, int64_t> cfem_dof_global2local_map;

  // Local cells first
  std::set<uint64_t> partition_bndry_vertex_id_set;
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns[u].num_components;
        for (size_t c = 0; c < num_components; ++c)
        {
          const int64_t dof_dfem_map = pwld_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          const int64_t dof_cfem_map = pwlc_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          const int64_t dof_cfem_map_global = pwlc_sdm.MapDOF(cell, i, uk_man, u, c);

          cfem_dof_global2local_map[dof_cfem_map_global] = dof_cfem_map;

          const double phi_value = input[dof_dfem_map];

          cont_input[dof_cfem_map] += phi_value;
          cont_input_ctr[dof_cfem_map] += 1.0;
        } // for component c
      }   // for unknown u
    }     // for node i

    for (const auto& face : cell.faces)
      if (face.has_neighbor)
        if (not grid->IsCellLocal(face.neighbor_id))
          for (const uint64_t vid : face.vertex_ids)
            partition_bndry_vertex_id_set.insert(vid);
  } // for local cell

  // Ghost cells
  const auto ghost_cell_ids = grid->cells.GetGhostGlobalIDs();
  const auto& vid_set = partition_bndry_vertex_id_set;
  for (const auto global_id : ghost_cell_ids)
  {
    const auto& cell = grid->cells[global_id];
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      if (vid_set.find(cell.vertex_ids[i]) == vid_set.end())
        continue;

      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns[u].num_components;
        for (size_t c = 0; c < num_components; ++c)
        {
          const int64_t dof_dfem_map_global = pwld_sdm.MapDOF(cell, i, uk_man, u, c);
          const int64_t dof_cfem_map_global = pwlc_sdm.MapDOF(cell, i, uk_man, u, c);
          if (cfem_dof_global2local_map.count(dof_cfem_map_global) > 0)
          {
            const int64_t dof_dfem_map = dfem_dof_global2local_map.at(dof_dfem_map_global);
            const int64_t dof_cfem_map = cfem_dof_global2local_map[dof_cfem_map_global];

            const double phi_value = input_with_ghosts[dof_dfem_map];

            cont_input[dof_cfem_map] += phi_value;
            cont_input_ctr[dof_cfem_map] += 1.0;
          }
        } // for component
      }   // for unknown
    }     // for node i
  }       // for ghost cell

  // Compute nodal averages
  {
    const size_t num_vals = cont_input.size();
    for (size_t k = 0; k < num_vals; ++k)
      cont_input[k] /= cont_input_ctr[k];
  }

  // Project back to dfem
  output = input;
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = pwld_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_components = uk_man.unknowns[u].num_components;
        for (size_t c = 0; c < num_components; ++c)
        {
          const int64_t dof_dfem_map = pwld_sdm.MapDOFLocal(cell, i, uk_man, u, c);
          const int64_t dof_cfem_map = pwlc_sdm.MapDOFLocal(cell, i, uk_man, u, c);

          const double phi_value = cont_input[dof_cfem_map];

          output[dof_dfem_map] = phi_value;
        } // for component c
      }   // for unknown u
    }     // for node i
  }       // for local cell
}

void
LBSKEigenAcceleration::CopyOnlyPhi0(const std::vector<double>& phi_in,
                                    std::vector<double>& phi_local)
{
  const auto& lbs_sdm = do_problem_.GetSpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();
  const auto& phi_uk_man = do_problem_.GetUnknownManager();
  const int gsi = front_gs_.groups.front().id;
  const size_t gss = front_gs_.groups.size();
  const size_t diff_num_local_dofs = pwlc_ptr_ ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                                               : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  if (pwlc_ptr_)
    NodallyAveragedPWLDVector(phi_in, copy_only_phi0_tmp_);
  const std::vector<double>* const phi_data = pwlc_ptr_ ? &copy_only_phi0_tmp_ : &phi_in;

  phi_local.resize(diff_num_local_dofs);
  std::fill(phi_local.begin(), phi_local.end(), 0.0);

  for (const auto& cell : do_problem_.GetGrid()->local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      const int64_t lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);
      const auto input_begin = phi_data->begin() + lbs_phi_map;
      const auto output_begin = phi_local.begin() + diff_phi_map;
      std::copy_n(input_begin, gss, output_begin);
    } // for node
  }   // for cell

  if (pwlc_ptr_)
    copy_only_phi0_tmp_.clear();
}

void
LBSKEigenAcceleration::ProjectBackPhi0(const std::vector<double>& input,
                                       std::vector<double>& output) const
{
  const auto& lbs_sdm = do_problem_.GetSpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();
  const auto& phi_uk_man = do_problem_.GetUnknownManager();
  const int gsi = front_gs_.groups.front().id;
  const size_t gss = front_gs_.groups.size();
  const size_t diff_num_local_dofs = pwlc_ptr_ ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                                               : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  OpenSnLogicalErrorIf(input.size() != diff_num_local_dofs, "Input vector size mismatch");

  const auto output_size = lbs_sdm.GetNumLocalDOFs(phi_uk_man);
  if (output.empty())
    output.resize(output_size, 0.0);
  OpenSnLogicalErrorIf(output.size() != output_size, "Output vector size mismatch");

  for (const auto& cell : do_problem_.GetGrid()->local_cells)
  {
    const auto& cell_mapping = lbs_sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t diff_phi_map = diff_sdm.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      const int64_t lbs_phi_map = lbs_sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);
      const auto input_begin = input.begin() + diff_phi_map;
      const auto output_begin = output.begin() + lbs_phi_map;
      std::copy_n(input_begin, gss, output_begin);
    } // for dof
  }   // for cell
}

} // namespace opensn
