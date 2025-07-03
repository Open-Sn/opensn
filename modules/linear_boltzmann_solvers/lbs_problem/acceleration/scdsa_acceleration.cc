// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/diffusion_pwlc_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/scdsa_acceleration.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/vector_ghost_communicator/vector_ghost_communicator.h"
#include "framework/object_factory.h"

namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, SCDSAAcceleration);

InputParameters
SCDSAAcceleration::GetInputParameters()
{
  auto params = LBSKEigenAcceleration::GetInputParameters();

  params.AddOptionalParameter("sdm", "pwld", "Spatial discretization");
  params.AddOptionalParameter(
    "pi_max_its", 50, "Maximum allowable iterations for inner power iterations");
  params.AddOptionalParameter(
    "pi_k_tol", 1.0e-10, "k-eigenvalue tolerance for the inner power iterations");

  params.ConstrainParameterRange("sdm", AllowableRangeList::New({"pwld", "pwlc"}));

  params.ChangeExistingParamToOptional("name", "SCDSAAcceleration");

  return params;
}

std::shared_ptr<SCDSAAcceleration>
SCDSAAcceleration::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SCDSAAcceleration>("lbs::SCDSAAcceleration", params);
}

SCDSAAcceleration::SCDSAAcceleration(const InputParameters& params)
  : LBSKEigenAcceleration(params),
    sdm_(params.GetParamValue<std::string>("sdm")),
    pi_max_its_(params.GetParamValue<int>("pi_max_its")),
    pi_k_tol_(params.GetParamValue<double>("pi_k_tol")),
    front_gs_(groupsets_.front())
{
  if (lbs_problem_.GetGroupsets().size() != 1)
    throw std::logic_error("The SCDSA acceleration scheme is only implemented for "
                           "problems with a single groupset.");

  // If using the AAH solver with one sweep, a few iterations need to be done
  // to get rid of the junk in the unconverged lagged angular fluxes.  Five
  // sweeps is a guess at how many initial sweeps are necessary.
  if (const auto do_problem = dynamic_cast<DiscreteOrdinatesProblem*>(&lbs_problem_))
    if (do_problem->GetSweepType() == "AAH" and front_gs_.max_iterations == 1)
      throw std::logic_error("The AAH solver is not stable for single-sweep methods due to "
                             "the presence of lagged angular fluxes.  Multiple sweeps are "
                             "allowed, however, the number of sweeps required to get sensible "
                             "results is not well studied and problem dependent.");
}

void
SCDSAAcceleration::Initialize()
{
  front_wgs_solver_ = lbs_problem_.GetWGSSolvers().at(front_gs_.id);
  front_wgs_context_ = std::dynamic_pointer_cast<WGSContext>(front_wgs_solver_->GetContext());
  OpenSnLogicalErrorIf(not front_wgs_context_, ": Casting failed");

  // Make UnknownManager
  const size_t num_gs_groups = front_gs_.groups.size();
  UnknownManager uk_man;
  uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

  // Make boundary conditions
  const auto bcs = TranslateBCs(lbs_problem_.GetSweepBoundaries(), true);

  // Make xs map
  const auto matid_2_mgxs_map = PackGroupsetXS(
    lbs_problem_.GetMatID2XSMap(), front_gs_.groups.front().id, front_gs_.groups.back().id);

  // Create solver
  const auto& sdm = lbs_problem_.GetSpatialDiscretization();
  const auto& unit_cell_matrices = lbs_problem_.GetUnitCellMatrices();

  if (sdm_ == "pwld")
  {
    diffusion_solver_ = std::make_shared<DiffusionMIPSolver>(std::string(GetName() + "_WGDSA"),
                                                             sdm,
                                                             uk_man,
                                                             bcs,
                                                             matid_2_mgxs_map,
                                                             unit_cell_matrices,
                                                             false,
                                                             true);
  }
  else
  {
    continuous_sdm_ptr_ = PieceWiseLinearContinuous::New(sdm.GetGrid());
    diffusion_solver_ = std::make_shared<DiffusionPWLCSolver>(std::string(GetName() + "_WGDSA"),
                                                              *continuous_sdm_ptr_,
                                                              uk_man,
                                                              bcs,
                                                              matid_2_mgxs_map,
                                                              unit_cell_matrices,
                                                              false,
                                                              true);
    lbs_pwld_ghost_info_ = MakePWLDVecGhostCommInfo();

    const auto& cfem_sdm = *continuous_sdm_ptr_;
    const auto ghost_dof_ids = cfem_sdm.GetGhostDOFIndices(lbs_problem_.GetUnknownManager());
  }

  auto& ds = *diffusion_solver_;
  ds.options.residual_tolerance = l_abs_tol_;
  ds.options.max_iters = max_iters_;
  ds.options.verbose = verbose_;
  ds.options.additional_options_string = petsc_options_;

  log.Log() << "Initializing diffusion solver";
  ds.Initialize();
  opensn::mpi_comm.barrier();
  log.Log() << "Done Initializing diffusion solver";

  log.Log() << "Assembling A and b";
  std::vector<double> dummy_rhs;
  if (sdm_ == "pwld")
    dummy_rhs.assign(sdm.GetNumLocalDOFs(uk_man), 0.0);
  else
    dummy_rhs.assign(continuous_sdm_ptr_->GetNumLocalAndGhostDOFs(uk_man), 0.0);
  ds.AssembleAand_b(dummy_rhs);
  log.Log() << "Done Assembling A and b";
}

void
SCDSAAcceleration::PreExecute()
{
  phi_temp_ = phi_old_local_;

  // If using multiple sweeps, the algorithm requires a delta-phi between two
  // successive iterations. This is not possible to obtain with the standard solve
  // routine.  Instead, a follow-up sweep must be performed.  If using two or more
  // sweeps, reduce the sweep count by one to ensure the user gets the requested
  // number of sweeps per outer.
  extra_sweep_ = front_gs_.max_iterations > 1;
  if (extra_sweep_)
    front_gs_.max_iterations -= 1;
}

void
SCDSAAcceleration::PrePowerIteration()
{
  Sf_ell_ = q_moments_local_;
  CopyOnlyPhi0(Sf_ell_, Sf0_ell_);

  // Init old scalar flux storage. If using a single sweep, set
  // it to the existing scalar flux
  if (!extra_sweep_)
    CopyOnlyPhi0(phi_old_local_, phi0_ell_);
}

double
SCDSAAcceleration::PostPowerIteration()
{
  // Lambda for the creation of scattering sources but the
  // input vector is only the zeroth moment
  auto SetLBSScatterSourcePhi0 = [this](const std::vector<double>& input,
                                        const bool additive,
                                        const bool suppress_wg_scat = false)
  {
    ProjectBackPhi0(input, phi_temp_);
    solver_->SetLBSScatterSource(phi_temp_, additive, suppress_wg_scat);
  };

  if (extra_sweep_)
  {
    // Set the old scalar flux to the scalar flux before the last sweep
    CopyOnlyPhi0(phi_new_local_, phi0_ell_);

    // Do an extra sweep
    q_moments_local_ = Sf_ell_;
    solver_->SetLBSScatterSource(phi_old_local_, true);
    front_wgs_context_->ApplyInverseTransportOperator(SourceFlags());
  }

  // Store the intermediate scalar flux
  CopyOnlyPhi0(phi_new_local_, phi0_star_);

  // Power Iteration Acceleration
  SetLBSScatterSourcePhi0(phi0_star_ - phi0_ell_, false);
  CopyOnlyPhi0(q_moments_local_, Ss0_res_);

  double production_k = lbs_problem_.ComputeFissionProduction(phi_new_local_);

  for (auto& vec : {&epsilon_k_, &epsilon_kp1_})
  {
    vec->resize(phi0_star_.size());
    std::fill(vec->begin(), vec->end(), 0.0);
  }

  double lambda_k = solver_->GetEigenvalue();
  double lambda_kp1 = lambda_k;

  for (size_t k = 0; k < pi_max_its_; ++k)
  {
    ProjectBackPhi0(epsilon_k_ + phi0_star_, phi_temp_);
    solver_->SetLBSFissionSource(phi_temp_, false);
    Scale(q_moments_local_, 1.0 / lambda_k);

    CopyOnlyPhi0(q_moments_local_, Sf0_aux_);

    // Inner iterations seem extremely wasteful. Set this to 1 iteration here for further
    // investigation.
    for (int i = 0; i < 1; ++i)
    {
      SetLBSScatterSourcePhi0(epsilon_k_, false, true);

      CopyOnlyPhi0(q_moments_local_, Ss0_);

      // Solve the diffusion system
      diffusion_solver_->Assemble_b(Ss0_ + Sf0_aux_ + Ss0_res_ - Sf0_ell_);
      diffusion_solver_->Solve(epsilon_kp1_, true);

      epsilon_k_ = epsilon_kp1_;
    }

    ProjectBackPhi0(epsilon_kp1_ + phi0_star_, phi_old_local_);

    const double production_kp1 = lbs_problem_.ComputeFissionProduction(phi_old_local_);

    lambda_kp1 = production_kp1 / (production_k / lambda_k);

    const double lambda_change = std::fabs(lambda_kp1 / lambda_k - 1.0);
    if (verbose_ >= 1)
      log.Log() << "PISCDSA iteration " << k << " lambda " << lambda_kp1 << " lambda change "
                << lambda_change;

    if (lambda_change < pi_k_tol_)
      break;

    lambda_k = lambda_kp1;
    epsilon_k_ = epsilon_kp1_;
    production_k = production_kp1;
  } // acceleration loop

  ProjectBackPhi0(epsilon_kp1_ + phi0_star_, phi_new_local_);
  LBSVecOps::GSScopedCopyPrimarySTLvectors(
    lbs_problem_, front_gs_, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);

  return lambda_kp1;
}

SCDSAAcceleration::GhostInfo
SCDSAAcceleration::MakePWLDVecGhostCommInfo() const
{
  const auto& sdm = lbs_problem_.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem_.GetUnknownManager();

  log.Log() << "Making PWLD ghost communicator";

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(uk_man);
  const size_t num_global_dofs = sdm.GetNumGlobalDOFs(uk_man);

  log.Log() << "Number of global dofs" << num_global_dofs;

  const size_t num_unknowns = uk_man.unknowns.size();

  // Build a list of global ids
  std::set<int64_t> global_dof_ids_set;

  const auto& grid = lbs_problem_.GetGrid();
  const auto ghost_cell_ids = grid->cells.GetGhostGlobalIDs();
  for (const auto global_id : ghost_cell_ids)
  {
    const auto& cell = grid->cells[global_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        const size_t num_comps = uk_man.unknowns[u].num_components;
        for (size_t c = 0; c < num_comps; ++c)
        {
          const int64_t dof_map = sdm.MapDOF(cell, i, uk_man, u, c);
          global_dof_ids_set.insert(dof_map);
        } // for component
      }   // for unknown
    }     // for node i
  }       // for ghost cell

  // Convert the list to a vector
  std::vector<int64_t> global_indices(global_dof_ids_set.begin(), global_dof_ids_set.end());

  // Create the vector ghost communicator
  auto vgc = std::make_shared<VectorGhostCommunicator>(
    num_local_dofs, num_global_dofs, global_indices, mpi_comm);

  // Create the map
  std::map<int64_t, int64_t> ghost_global_id_2_local_map;
  int64_t k = 0;
  for (const auto ghost_id : global_indices)
  {
    ghost_global_id_2_local_map[ghost_id] = static_cast<int64_t>(num_local_dofs + k++);
  }

  log.Log() << "Done making PWLD ghost communicator";
  return {vgc, ghost_global_id_2_local_map};
}

void
SCDSAAcceleration::CopyOnlyPhi0(const std::vector<double>& phi_in, std::vector<double>& phi_local)
{
  const auto& lbs_sdm = lbs_problem_.GetSpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();
  const auto& phi_uk_man = lbs_problem_.GetUnknownManager();
  const int gsi = front_gs_.groups.front().id;
  const size_t gss = front_gs_.groups.size();
  const size_t diff_num_local_dofs = continuous_sdm_ptr_
                                       ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                                       : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  if (continuous_sdm_ptr_)
    NodallyAveragedPWLDVector(phi_in, copy_only_phi0_tmp_);
  const std::vector<double>* const phi_data = continuous_sdm_ptr_ ? &copy_only_phi0_tmp_ : &phi_in;

  phi_local.resize(diff_num_local_dofs);
  std::fill(phi_local.begin(), phi_local.end(), 0.0);

  for (const auto& cell : lbs_problem_.GetGrid()->local_cells)
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

  if (continuous_sdm_ptr_)
    copy_only_phi0_tmp_.clear();
}

void
SCDSAAcceleration::ProjectBackPhi0(const std::vector<double>& input,
                                   std::vector<double>& output) const
{
  const auto& lbs_sdm = lbs_problem_.GetSpatialDiscretization();
  const auto& diff_sdm = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();
  const auto& phi_uk_man = lbs_problem_.GetUnknownManager();
  const int gsi = front_gs_.groups.front().id;
  const size_t gss = front_gs_.groups.size();
  const size_t diff_num_local_dofs = continuous_sdm_ptr_
                                       ? diff_sdm.GetNumLocalAndGhostDOFs(diff_uk_man)
                                       : diff_sdm.GetNumLocalDOFs(diff_uk_man);

  OpenSnLogicalErrorIf(input.size() != diff_num_local_dofs, "Vector size mismatch");

  for (const auto& cell : lbs_problem_.GetGrid()->local_cells)
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

void
SCDSAAcceleration::NodallyAveragedPWLDVector(const std::vector<double>& input,
                                             std::vector<double>& output) const
{
  const auto& pwld_sdm = lbs_problem_.GetSpatialDiscretization();
  const auto& pwlc_sdm = diffusion_solver_->GetSpatialDiscretization();
  const auto& uk_man = lbs_problem_.GetUnknownManager();

  const auto& vgc = lbs_pwld_ghost_info_.vector_ghost_communicator;
  const auto& dfem_dof_global2local_map = lbs_pwld_ghost_info_.ghost_global_id_2_local_map;

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
} // namespace opensn
