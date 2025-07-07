// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/diffusion_mip_solver.h"
#include "caliper/cali.h"

namespace opensn
{

void
WGDSA::Init(LBSProblem& lbs_problem, LBSGroupset& groupset, bool vaccum_bcs_are_dirichlet)
{
  CALI_CXX_MARK_SCOPE("WGDSA::Init");

  if (groupset.apply_wgdsa)
  {
    // Make UnknownManager
    const size_t num_gs_groups = groupset.groups.size();
    opensn::UnknownManager uk_man;
    uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

    // Make boundary conditions
    auto sweep_boundaries = lbs_problem.GetSweepBoundaries();
    auto bcs = TranslateBCs(sweep_boundaries, vaccum_bcs_are_dirichlet);

    // Make xs map
    auto block_id_to_xs_map = lbs_problem.GetMatID2XSMap();
    auto matid_2_mgxs_map =
      PackGroupsetXS(block_id_to_xs_map, groupset.groups.front().id, groupset.groups.back().id);

    // Create solver
    const auto& sdm = lbs_problem.GetSpatialDiscretization();
    auto lbs_name = lbs_problem.GetName();
    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(lbs_name + "_WGDSA"),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       lbs_problem.GetUnitCellMatrices(),
                                                       false,
                                                       true);
    ParameterBlock block;

    solver->options.residual_tolerance = groupset.wgdsa_tol;
    solver->options.max_iters = groupset.wgdsa_max_iters;
    solver->options.verbose = groupset.wgdsa_verbose;
    solver->options.additional_options_string = groupset.wgdsa_string;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.wgdsa_solver = solver;
  }
}

void
WGDSA::AssembleDeltaPhiVector(LBSProblem& lbs_problem,
                              const LBSGroupset& groupset,
                              const std::vector<double>& phi_in,
                              std::vector<double>& delta_phi_local)
{
  CALI_CXX_MARK_SCOPE("WGDSA::AssembleDeltaPhiVector");

  const auto grid = lbs_problem.GetGrid();
  const auto& sdm = lbs_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = lbs_problem.GetUnknownManager();
  const auto& block_id_to_xs_map = lbs_problem.GetMatID2XSMap();

  const int gsi = groupset.groups.front().id;
  const size_t gss = groupset.groups.size();

  delta_phi_local.clear();
  delta_phi_local.assign(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto& sigma_s = block_id_to_xs_map.at(cell.block_id)->GetSigmaSGtoG();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* delta_phi_mapped = &delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; ++g)
      {
        delta_phi_mapped[g] = sigma_s[gsi + g] * phi_in_mapped[g];
      }
    }
  }
}

void
WGDSA::DisassembleDeltaPhiVector(LBSProblem& lbs_problem,
                                 const LBSGroupset& groupset,
                                 const std::vector<double>& delta_phi_local,
                                 std::vector<double>& ref_phi_new)
{
  CALI_CXX_MARK_SCOPE("WGDSA::DisassembleDeltaPhiVector");

  const auto grid = lbs_problem.GetGrid();
  const auto& sdm = lbs_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = lbs_problem.GetUnknownManager();

  const int gsi = groupset.groups.front().id;
  const size_t gss = groupset.groups.size();

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* delta_phi_mapped = &delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (int g = 0; g < gss; ++g)
        phi_new_mapped[g] += delta_phi_mapped[g];
    }
  }
}

void
WGDSA::CleanUp(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("WGDSA::CleanUp");

  if (groupset.apply_wgdsa)
    groupset.wgdsa_solver = nullptr;
}

std::vector<double>
WGDSA::WGSCopyOnlyPhi0(LBSProblem& lbs_problem,
                       const LBSGroupset& groupset,
                       const std::vector<double>& phi_in)
{
  CALI_CXX_MARK_SCOPE("WGDSA::WGSCopyOnlyPhi0");

  const auto grid = lbs_problem.GetGrid();
  const auto& sdm = lbs_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = lbs_problem.GetUnknownManager();

  const int gsi = groupset.groups.front().id;
  const size_t gss = groupset.groups.size();

  std::vector<double> output_phi_local(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; ++g)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell

  return output_phi_local;
}

void
WGDSA::GSProjectBackPhi0(LBSProblem& lbs_problem,
                         const LBSGroupset& groupset,
                         const std::vector<double>& input,
                         std::vector<double>& output)
{
  CALI_CXX_MARK_SCOPE("WGDSA::GSProjectBackPhi0");

  const auto grid = lbs_problem.GetGrid();
  const auto& sdm = lbs_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = lbs_problem.GetUnknownManager();

  const int gsi = groupset.groups.front().id;
  const size_t gss = groupset.groups.size();

  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[dphi_map];
      double* output_mapped = &output[phi_map];

      for (int g = 0; g < gss; ++g)
        output_mapped[g] = input_mapped[g];
    } // for dof
  }   // for cell
}

} // namespace opensn
