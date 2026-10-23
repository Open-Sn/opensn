// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/tgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

void
TGDSA::Init(DiscreteOrdinatesProblem& do_problem, LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa)
  {
    CALI_CXX_MARK_SCOPE("Acceleration/TGDSA");

    // DiscreteOrdinatesProblem::CheckConfigurationErrors rejects single-group groupsets.
    OpenSnLogicalErrorIf(groupset.GetNumGroups() < 2,
                         do_problem.GetName() + ": TGDSA on single-group groupset " +
                           std::to_string(groupset.id) + ".");

    const auto& sdm = do_problem.GetSpatialDiscretization();
    const auto& uk_man = sdm.UNITARY_UNKNOWN_MANAGER;
    const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();
    const auto& sweep_boundaries = do_problem.GetSweepBoundaries();

    // Make boundary conditions
    auto bcs = TranslateBCs(sweep_boundaries);

    // Make TwoGridInfo
    groupset.tg_acceleration_info_.map_mat_id_2_tginfo.clear();
    for (const auto& mat_id_xs_pair : block_id_to_xs_map)
    {
      const auto& mat_id = mat_id_xs_pair.first;
      const auto& xs = mat_id_xs_pair.second;

      // The spectrum is that of Jacobi iteration with converged within-group scattering (JFULL),
      // with or without WGDSA. Without WGDSA, the spectrum of the single-sweep iteration
      // (JPARTIAL) matches the iteration more closely in a homogeneous medium, but it is
      // dominated by each material's most strongly self-scattering group, so neighboring
      // materials can have spectra in different groups. The collapsed correction is then mapped
      // to the wrong groups across material interfaces, and the iteration can diverge.
      TwoGridCollapsedInfo tginfo = MakeTwoGridCollapsedInfo(*xs,
                                                             EnergyCollapseScheme::JFULL,
                                                             groupset.first_group,
                                                             groupset.last_group,
                                                             do_problem.GetDSATimeAbsorptionScale(),
                                                             do_problem.GetOptions().adjoint);

      groupset.tg_acceleration_info_.map_mat_id_2_tginfo.insert(
        std::make_pair(mat_id, std::move(tginfo)));
    }

    // Make xs map
    std::map<unsigned int, Multigroup_D_and_sigR> matid_2_mgxs_map;
    for (const auto& matid_xs_pair : block_id_to_xs_map)
    {
      const auto& mat_id = matid_xs_pair.first;

      const auto& tg_info = groupset.tg_acceleration_info_.map_mat_id_2_tginfo.at(mat_id);

      matid_2_mgxs_map.insert(std::make_pair(
        mat_id, Multigroup_D_and_sigR{{tg_info.collapsed_D}, {tg_info.collapsed_sig_a}}));
    }

    // Create solver
    const auto lbs_name = do_problem.GetName();
    const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(lbs_name + "_TGDSA"),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       unit_cell_matrices,
                                                       false,
                                                       true);

    solver->options.residual_tolerance = groupset.tgdsa_tol;
    solver->options.max_iters = groupset.tgdsa_max_iters;
    solver->options.solver_policy = groupset.tgdsa_solver_policy;
    solver->options.direct_solve_threshold = groupset.tgdsa_direct_solve_threshold;
    solver->options.verbose = groupset.tgdsa_verbose;
    solver->options.additional_options_string = groupset.tgdsa_string;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.tgdsa_solver = solver;
  }
}

void
TGDSA::AssembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                              const LBSGroupset& groupset,
                              const std::vector<double>& phi_in,
                              std::vector<double>& delta_phi_local)
{

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& phi_uk_man = do_problem.GetUnknownManager();
  const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();
  // The residual is the scattering between groups of the groupset, matching the JFULL spectrum
  // (see TGDSA::Init). Groups outside the groupset are not part of its iteration error.

  auto local_node_count = do_problem.GetLocalNodeCount();
  if (delta_phi_local.size() != local_node_count)
    delta_phi_local.assign(local_node_count, 0.0);
  else
    std::fill(delta_phi_local.begin(), delta_phi_local.end(), 0.0);

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto& S =
      block_id_to_xs_map.at(cell->block_id)->GetTransferMatrix(0, do_problem.GetOptions().adjoint);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i);
      const auto phi_map = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, 0);

      double& delta_phi_mapped = delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
      {
        double R_g = 0.0;
        for (const auto& [row_g, gprime, sigma_sm] : S.Row(gsi + g))
          if (gprime >= gsi and gprime < gsi + gss and gprime != (gsi + g))
            R_g += sigma_sm * phi_in_mapped[gprime];

        delta_phi_mapped += R_g;
      }
    }
  }
}

void
TGDSA::DisassembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                                 const LBSGroupset& groupset,
                                 const std::vector<double>& delta_phi_local,
                                 std::vector<double>& ref_phi_new)
{

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& phi_uk_man = do_problem.GetUnknownManager();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  const auto& map_mat_id_2_tginfo = groupset.tg_acceleration_info_.map_mat_id_2_tginfo;

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    const auto& xi_g = map_mat_id_2_tginfo.at(cell->block_id).spectrum;

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i);
      const auto phi_map = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi);

      const double delta_phi_mapped = delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
        phi_new_mapped[g] += delta_phi_mapped * xi_g[gsi + g];
    }
  }
}

void
TGDSA::CleanUp(LBSGroupset& groupset)
{

  if (groupset.apply_tgdsa)
  {
    groupset.tgdsa_solver = nullptr;
    groupset.tg_acceleration_info_.map_mat_id_2_tginfo.clear();
  }
}

} // namespace opensn
