// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/tgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "caliper/cali.h"

namespace opensn
{

void
TGDSA::Init(DiscreteOrdinatesProblem& do_problem, LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("TGDSA::Init");

  if (groupset.apply_tgdsa)
  {
    const auto& sdm = do_problem.GetSpatialDiscretization();
    const auto& uk_man = sdm.UNITARY_UNKNOWN_MANAGER;
    const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();
    const auto& sweep_boundaries = do_problem.GetSweepBoundaries();

    // Make boundary conditions
    auto bcs = TranslateBCs(sweep_boundaries);

    // Make TwoGridInfo
    for (const auto& mat_id_xs_pair : block_id_to_xs_map)
    {
      const auto& mat_id = mat_id_xs_pair.first;
      const auto& xs = mat_id_xs_pair.second;

      TwoGridCollapsedInfo tginfo = MakeTwoGridCollapsedInfo(*xs, EnergyCollapseScheme::JFULL);

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
  CALI_CXX_MARK_SCOPE("TGDSA::AssembleDeltaPhiVector");

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& phi_uk_man = do_problem.GetUnknownManager();
  const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  auto local_node_count = do_problem.GetLocalNodeCount();
  delta_phi_local.clear();
  delta_phi_local.assign(local_node_count, 0.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto& S = block_id_to_xs_map.at(cell.block_id)->GetTransferMatrix(0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(cell, i);
      const auto phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

      double& delta_phi_mapped = delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
      {
        double R_g = 0.0;
        for (const auto& [row_g, gprime, sigma_sm] : S.Row(gsi + g))
          if (gprime >= gsi and gprime != (gsi + g))
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
  CALI_CXX_MARK_SCOPE("TGDSA::DisassembleDeltaPhiVector");

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& phi_uk_man = do_problem.GetUnknownManager();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  const auto& map_mat_id_2_tginfo = groupset.tg_acceleration_info_.map_mat_id_2_tginfo;

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    const auto& xi_g = map_mat_id_2_tginfo.at(cell.block_id).spectrum;

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(cell, i);
      const auto phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

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
  CALI_CXX_MARK_SCOPE("TGDSA::CleanUp");

  if (groupset.apply_tgdsa)
    groupset.tgdsa_solver = nullptr;
}

} // namespace opensn
