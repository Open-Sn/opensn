// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/csda_sweep_helper.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "caliper/cali.h"
#include <algorithm>
#include <string>

namespace opensn
{

namespace
{

void
SetSolverOptions(DiffusionMIPSolver& solver, const LBSGroupset& groupset)
{
  solver.options.residual_tolerance = groupset.wgdsa_tol;
  solver.options.max_iters = groupset.wgdsa_max_iters;
  solver.options.solver_policy = groupset.wgdsa_solver_policy;
  solver.options.direct_solve_threshold = groupset.wgdsa_direct_solve_threshold;
  solver.options.verbose = groupset.wgdsa_verbose;
  solver.options.additional_options_string = groupset.wgdsa_string;
  solver.options.h_perp_convention = groupset.wgdsa_h_perp_convention;
  if (DSABoundaryClosureFromString(groupset.wgdsa_boundary_closure) ==
      DSABoundaryClosure::CLOSURE_R)
    solver.SetBoundaryClosureQuadrature(groupset.quadrature->GetOmegas(),
                                        groupset.quadrature->GetWeights());
}

/// True when a reflecting boundary of the groupset lags its incoming angular flux by one sweep.
bool
HasDelayedReflectingBoundary(const DiscreteOrdinatesProblem& do_problem,
                             const LBSGroupset& groupset)
{
  const auto& boundaries = do_problem.GetSweepBoundaries();
  return std::any_of(boundaries.begin(),
                     boundaries.end(),
                     [&groupset](const auto& bid_bndry)
                     {
                       const auto& bndry = *bid_bndry.second;
                       return bndry.IsReflecting() and
                              bndry.CountDelayedAngularDOFsNew(groupset.id) > 0;
                     });
}

/**
 * Passes to `add(cell, i, g, value)` the boundary source of the WGDSA correction on the reflecting
 * boundaries that lag their incoming angular flux. There the incoming error is the reflected
 * outgoing error plus the change delta psi of the lagged flux that the sweep did not see, so the
 * correction satisfies D n.grad(delta phi) = J, with J the inward partial current of delta psi
 * (SweepBoundary::TallyDelayedIncomingCurrent). With J interpolated by the face basis, the source
 * of node i is the face term sum_j (int_f b_i b_j) J_j.
 */
template <typename AddFn>
void
TallyDelayedReflectingSource(const DiscreteOrdinatesProblem& do_problem,
                             const LBSGroupset& groupset,
                             const bool subtract_old,
                             const AddFn& add)
{
  const auto& grid = *do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();

  const auto tally = [&](std::uint32_t cell_local_id,
                         unsigned int face,
                         unsigned int face_node,
                         unsigned int group_idx,
                         double current)
  {
    if (current == 0.0)
      return;
    const auto& cell = *grid.GetLocalCells()[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& face_mass = unit_cell_matrices[cell_local_id].intS_shapeI_shapeJ[face];
    const auto j = cell_mapping.MapFaceNode(face, face_node);
    for (size_t fi = 0; fi < cell_mapping.GetNumFaceNodes(face); ++fi)
    {
      const auto i = cell_mapping.MapFaceNode(face, fi);
      add(cell, i, group_idx, face_mass(i, j) * current);
    }
  };

  for (const auto& [bid, bndry] : do_problem.GetSweepBoundaries())
    bndry->TallyDelayedIncomingCurrent(groupset.id, *groupset.quadrature, subtract_old, tally);
}

/// True when CSDA is enabled and the groupset contains a charged-particle group.
bool
UsesCSDAOperator(const DiscreteOrdinatesProblem& do_problem, const LBSGroupset& groupset)
{
  if (not do_problem.GetOptions().csda_enabled)
    return false;
  const auto ranges =
    FindCSDAProblemChargedGroupRanges(do_problem.GetBlockID2XSMap(), do_problem.GetNumGroups());
  return std::any_of(
    ranges.begin(),
    ranges.end(),
    [&groupset](const auto& range)
    { return range.first <= groupset.last_group and range.second > groupset.first_group; });
}

/**
 * Builds one single-group diffusion solver per groupset group. Charged groups add the CSDA
 * removal S_g/dE_g to the removal cross section and to the transport cross section of the
 * diffusion coefficient.
 */
void
InitCSDA(DiscreteOrdinatesProblem& do_problem,
         LBSGroupset& groupset,
         const std::map<uint64_t, BoundaryCondition>& bcs)
{
  const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();
  const auto energy =
    MultiGroupXS::ResolveEnergyGroupStructure(block_id_to_xs_map, do_problem.GetNumGroups());

  opensn::UnknownManager uk_man;
  uk_man.AddUnknown(UnknownType::VECTOR_N, 1);

  const auto& sdm = do_problem.GetSpatialDiscretization();
  std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);
  for (auto g = groupset.first_group; g <= groupset.last_group; ++g)
  {
    auto matid_2_mgxs_map =
      PackGroupsetXS(block_id_to_xs_map, g, g, do_problem.GetDSATimeAbsorptionScale());
    for (auto& [block_id, xs_g] : matid_2_mgxs_map)
    {
      const auto csda_data = MakeCSDAMaterialData(*block_id_to_xs_map.at(block_id), energy);
      if (not IsCSDAActive(csda_data, g))
        continue;
      const double removal = (*csda_data.stopping_power)[g] / (*csda_data.delta_e)[g];
      xs_g.Dg[0] = 1.0 / (1.0 / xs_g.Dg[0] + 3.0 * removal);
      xs_g.sigR[0] += removal;
    }

    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(do_problem.GetName()) +
                                                         "_WGDSA_g" + std::to_string(g),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       do_problem.GetUnitCellMatrices(),
                                                       false,
                                                       true);
    SetSolverOptions(*solver, groupset);
    solver->Initialize();
    solver->AssembleAand_b(dummy_rhs);
    groupset.wgdsa_csda_solvers.push_back(solver);
  }
}

/**
 * CSDA counterpart of the assemble-solve-disassemble sequence. The residual source of group g is
 * the P0 scattering of `phi_change` from every groupset group into g. The groups are then solved in
 * energy order; group g also receives the in-groupset scattering and CSDA slowing-down of the
 * corrections of groups above it, so the correction is coupled in energy like the sweep.
 */
void
ApplyCSDA(DiscreteOrdinatesProblem& do_problem,
          const LBSGroupset& groupset,
          const std::vector<double>& phi_change,
          const bool subtract_old_delayed_psi,
          std::vector<double>& phi_new)
{
  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& phi_uk_man = do_problem.GetUnknownManager();
  const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();
  const auto energy =
    MultiGroupXS::ResolveEnergyGroupStructure(block_id_to_xs_map, do_problem.GetNumGroups());
  const auto& solvers = groupset.wgdsa_csda_solvers;
  const auto& dphi_uk_man = solvers.front()->GetUnknownStructure();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();
  const auto num_dofs = sdm.GetNumLocalDOFs(dphi_uk_man);

  // Adds the P0 transfer into groupset group g from the groupset groups [0, g_end) of `source`.
  const auto add_transfer = [&](const SparseMatrix& S0,
                                const unsigned int g,
                                const unsigned int g_end,
                                const auto& source_at,
                                double& value)
  {
    const auto& cols = S0.rowI_indices[gsi + g];
    const auto& vals = S0.rowI_values[gsi + g];
    for (size_t t = 0; t < cols.size(); ++t)
      if (cols[t] >= gsi and cols[t] < gsi + g_end)
        value += vals[t] * source_at(cols[t] - gsi);
  };

  // Scattering residual of every groupset group. Each entry is replaced by its group's correction
  // once that group is solved.
  std::vector<std::vector<double>> correction(gss, std::vector<double>(num_dofs, 0.0));
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& transfer_matrices = block_id_to_xs_map.at(cell->block_id)->GetTransferMatrices();
    if (transfer_matrices.empty())
      continue;
    const size_t num_nodes = sdm.GetCellMapping(*cell).GetNumNodes();
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
      const double* dphi = &phi_change[sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi)];
      for (unsigned int g = 0; g < gss; ++g)
        add_transfer(
          transfer_matrices[0],
          g,
          gss,
          [dphi](unsigned int gp) { return dphi[gp]; },
          correction[g][dphi_map]);
    }
  }

  // Face source of each group on the reflecting boundaries that lag their incoming flux
  std::vector<std::vector<double>> boundary_source;
  if (HasDelayedReflectingBoundary(do_problem, groupset))
  {
    boundary_source.assign(gss, std::vector<double>(num_dofs, 0.0));
    TallyDelayedReflectingSource(
      do_problem,
      groupset,
      subtract_old_delayed_psi,
      [&](const Cell& cell, size_t i, unsigned int g, double value)
      { boundary_source[g][sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0)] += value; });
  }

  for (unsigned int g = 0; g < gss; ++g)
  {
    auto& q = correction[g];
    if (g > 0)
    {
      for (const auto& cell : grid->GetLocalCells())
      {
        const auto& xs = *block_id_to_xs_map.at(cell->block_id);
        const auto csda_data = MakeCSDAMaterialData(xs, energy);
        const bool slows_down =
          IsCSDAActive(csda_data, gsi + g) and IsCSDAActive(csda_data, gsi + g - 1);
        const double feed =
          slows_down ? (*csda_data.stopping_power)[gsi + g - 1] / (*csda_data.delta_e)[gsi + g - 1]
                     : 0.0;
        const auto& transfer_matrices = xs.GetTransferMatrices();
        const size_t num_nodes = sdm.GetCellMapping(*cell).GetNumNodes();
        for (size_t i = 0; i < num_nodes; ++i)
        {
          const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
          if (not transfer_matrices.empty())
            add_transfer(
              transfer_matrices[0],
              g,
              g,
              [&correction, dphi_map](unsigned int gp) { return correction[gp][dphi_map]; },
              q[dphi_map]);
          q[dphi_map] += feed * correction[g - 1][dphi_map];
        }
      }
    }
    solvers[g]->Assemble_b(q);
    if (not boundary_source.empty())
      solvers[g]->AddToRHS(boundary_source[g]);
    solvers[g]->Solve(q);
  }

  for (const auto& cell : grid->GetLocalCells())
  {
    const size_t num_nodes = sdm.GetCellMapping(*cell).GetNumNodes();
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
      double* phi_new_mapped = &phi_new[sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi)];
      for (unsigned int g = 0; g < gss; ++g)
        phi_new_mapped[g] += correction[g][dphi_map];
    }
  }
}

} // namespace

void
WGDSA::Init(DiscreteOrdinatesProblem& do_problem, LBSGroupset& groupset)
{
  if (groupset.apply_wgdsa)
  {
    CALI_CXX_MARK_SCOPE("Acceleration/WGDSA");

    // Make UnknownManager
    const auto num_gs_groups = groupset.GetNumGroups();
    opensn::UnknownManager uk_man;
    uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

    // Make boundary conditions
    const auto& sweep_boundaries = do_problem.GetSweepBoundaries();
    const auto closure = DSABoundaryClosureFromString(groupset.wgdsa_boundary_closure);
    auto bcs = TranslateBCs(sweep_boundaries, closure);

    if (UsesCSDAOperator(do_problem, groupset))
    {
      InitCSDA(do_problem, groupset, bcs);
      return;
    }

    // Make xs map
    const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();
    auto matid_2_mgxs_map = PackGroupsetXS(block_id_to_xs_map,
                                           groupset.first_group,
                                           groupset.last_group,
                                           do_problem.GetDSATimeAbsorptionScale());

    // Create solver
    const auto& sdm = do_problem.GetSpatialDiscretization();
    auto lbs_name = do_problem.GetName();
    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(lbs_name + "_WGDSA"),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       do_problem.GetUnitCellMatrices(),
                                                       false,
                                                       true);
    SetSolverOptions(*solver, groupset);

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.wgdsa_solver = solver;
  }
}

void
WGDSA::AssembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                              const LBSGroupset& groupset,
                              const std::vector<double>& phi_in,
                              std::vector<double>& delta_phi_local)
{

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = do_problem.GetUnknownManager();
  const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  const auto num_delta_phi_dofs = sdm.GetNumLocalDOFs(dphi_uk_man);
  if (delta_phi_local.size() != num_delta_phi_dofs)
    delta_phi_local.assign(num_delta_phi_dofs, 0.0);
  else
    std::fill(delta_phi_local.begin(), delta_phi_local.end(), 0.0);

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto& sigma_s = block_id_to_xs_map.at(cell->block_id)->GetSigmaSGtoG();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
      const auto phi_map = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi);

      double* delta_phi_mapped = &delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
      {
        delta_phi_mapped[g] = sigma_s[gsi + g] * phi_in_mapped[g];
      }
    }
  }
}

void
WGDSA::DisassembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                                 const LBSGroupset& groupset,
                                 const std::vector<double>& delta_phi_local,
                                 std::vector<double>& ref_phi_new)
{

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = do_problem.GetUnknownManager();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
      const auto phi_map = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi);

      const double* delta_phi_mapped = &delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
        phi_new_mapped[g] += delta_phi_mapped[g];
    }
  }
}

void
WGDSA::Apply(DiscreteOrdinatesProblem& do_problem,
             const LBSGroupset& groupset,
             const std::vector<double>& phi_change,
             const bool subtract_old_delayed_psi,
             std::vector<double>& phi_new)
{
  if (not groupset.wgdsa_csda_solvers.empty())
  {
    ApplyCSDA(do_problem, groupset, phi_change, subtract_old_delayed_psi, phi_new);
    return;
  }

  std::vector<double> delta_phi_local;
  AssembleDeltaPhiVector(do_problem, groupset, phi_change, delta_phi_local);
  groupset.wgdsa_solver->Assemble_b(delta_phi_local);
  if (HasDelayedReflectingBoundary(do_problem, groupset))
  {
    const auto& sdm = do_problem.GetSpatialDiscretization();
    const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
    std::vector<double> boundary_source(delta_phi_local.size(), 0.0);
    TallyDelayedReflectingSource(
      do_problem,
      groupset,
      subtract_old_delayed_psi,
      [&](const Cell& cell, size_t i, unsigned int g, double value)
      { boundary_source[sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, g)] += value; });
    groupset.wgdsa_solver->AddToRHS(boundary_source);
  }
  groupset.wgdsa_solver->Solve(delta_phi_local);
  DisassembleDeltaPhiVector(do_problem, groupset, delta_phi_local, phi_new);
}

void
WGDSA::CleanUp(LBSGroupset& groupset)
{

  if (groupset.apply_wgdsa)
  {
    groupset.wgdsa_solver = nullptr;
    groupset.wgdsa_csda_solvers.clear();
  }
}

std::vector<double>
WGDSA::WGSCopyOnlyPhi0(DiscreteOrdinatesProblem& do_problem,
                       const LBSGroupset& groupset,
                       const std::vector<double>& phi_in)
{

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = do_problem.GetUnknownManager();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  std::vector<double> output_phi_local(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
      const auto phi_map = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    } // for node
  } // for cell

  return output_phi_local;
}

void
WGDSA::GSProjectBackPhi0(DiscreteOrdinatesProblem& do_problem,
                         const LBSGroupset& groupset,
                         const std::vector<double>& input,
                         std::vector<double>& output)
{

  const auto grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& dphi_uk_man = groupset.wgdsa_solver->GetUnknownStructure();
  const auto& phi_uk_man = do_problem.GetUnknownManager();

  const auto gsi = groupset.first_group;
  const auto gss = groupset.GetNumGroups();

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto dphi_map = sdm.MapDOFLocal(*cell, i, dphi_uk_man, 0, 0);
      const auto phi_map = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[dphi_map];
      double* output_mapped = &output[phi_map];

      for (unsigned int g = 0; g < gss; ++g)
        output_mapped[g] = input_mapped[g];
    } // for dof
  } // for cell
}

} // namespace opensn
