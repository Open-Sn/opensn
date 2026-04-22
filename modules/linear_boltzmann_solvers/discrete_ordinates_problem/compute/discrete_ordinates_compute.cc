// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/utils/caliper_scopes.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <utility>
#include <unordered_set>

namespace opensn
{
namespace
{

std::vector<std::pair<unsigned int, unsigned int>>
FindChargedGroupRangesForBalance(const std::vector<double>& stopping_power)
{
  constexpr double tol = 1.0e-12;
  std::vector<std::pair<unsigned int, unsigned int>> ranges;

  unsigned int g = 0;
  while (g < stopping_power.size())
  {
    while (g < stopping_power.size() and std::abs(stopping_power[g]) <= tol)
      ++g;
    if (g >= stopping_power.size())
      break;

    const unsigned int g_begin = g;
    while (g < stopping_power.size() and std::abs(stopping_power[g]) > tol)
      ++g;
    ranges.emplace_back(g_begin, g);
  }

  return ranges;
}

} // namespace

namespace
{

double
GetInflow(const Cell& cell,
          const CellMapping& cell_mapping,
          const std::vector<UnitCellMatrices>& unit_cell_matrices,
          const CellFace& face,
          const std::shared_ptr<AngularQuadrature> quadrature,
          SweepBoundary& boundary,
          unsigned int f,
          int gs_id,
          unsigned int g)
{
  const auto& fe_vals = unit_cell_matrices.at(cell.local_id);
  const auto& int_f_shape_i = fe_vals.intS_shapeI[f];
  const unsigned int num_face_nodes = cell_mapping.GetNumFaceNodes(f);
  const size_t num_angles = quadrature->GetNumAngles();

  double inflow = 0.0;
  for (size_t n = 0; n < num_angles; ++n)
  {
    const auto& omega = quadrature->GetOmega(n);
    const auto& weight = quadrature->GetWeight(n);
    const double mu = omega.Dot(face.normal);

    if (mu >= 0.0)
      continue;

    for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
    {
      const unsigned int i = cell_mapping.MapFaceNode(f, fi);
      const double* psi_ptr = boundary.PsiIncoming(cell.local_id, f, fi, n, gs_id, g);
      const double psi_inc = (psi_ptr != nullptr) ? *psi_ptr : 0.0;
      inflow += mu * weight * int_f_shape_i(i) * psi_inc;
    }
  }

  return inflow;
}

} // namespace

BalanceTable
ComputeBalanceTable(DiscreteOrdinatesProblem& do_problem, double scaling_factor)
{
  opensn::mpi_comm.barrier();

  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  auto& phi_new_local = do_problem.GetPhiNewLocal();
  const auto& groupsets = do_problem.GetGroupsets();
  auto& q_moments_local = do_problem.GetQMomentsLocal();
  const auto saved_q_moments_local = q_moments_local;
  auto active_set_source_fn = do_problem.GetActiveSetSourceFunction();
  const auto& cell_transport_views = do_problem.GetCellTransportViews();
  const auto& cell_outflow_views = do_problem.GetCellOutflowViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& sweep_boundaries = do_problem.GetSweepBoundaries();
  const auto& block_id_to_xs_map = do_problem.GetBlockID2XSMap();
  const auto num_groups = do_problem.GetNumGroups();
  const auto& options = do_problem.GetOptions();
  const auto& phi_e_new_local = do_problem.GetPhiENewLocal();
  const auto time_dependent = do_problem.IsTimeDependent();
  const auto dt = time_dependent ? do_problem.GetTimeStep() : 0.0;
  const auto& psi_new_local = do_problem.GetPsiNewLocal();
  const auto& psi_old_local = do_problem.GetPsiOldLocal();
  const auto& first_collision_source = do_problem.GetFirstCollisionSourceMoments();

  // Get material source
  // This is done using the SetSource routine because it allows a lot of flexibility.
  auto mat_src = phi_new_local;
  mat_src.assign(mat_src.size(), 0.0);
  for (const auto& groupset : groupsets)
  {
    do_problem.ZeroQMoments();
    active_set_source_fn(groupset,
                         q_moments_local,
                         phi_new_local,
                         APPLY_FIXED_SOURCES | APPLY_AGS_FISSION_SOURCES |
                           APPLY_WGS_FISSION_SOURCES | APPLY_PREVIOUS_PRECURSOR_SOURCES);
    LBSVecOps::GSScopedCopyPrimarySTLvectors( // NOLINT(readability-suspicious-call-argument)
      do_problem,
      groupset,
      q_moments_local,
      mat_src);
  }

  // Compute absorption, material-source and in-flow
  double local_out_flow = 0.0;
  double local_in_flow = 0.0;
  double local_absorption = 0.0;
  double local_production = 0.0;
  double local_initial = 0.0;
  double local_final = 0.0;
  double local_csda_charge_deposition = 0.0;
  double local_csda_energy_deposition = 0.0;
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = discretization.GetCellMapping(*cell);
    const auto& transport_view = cell_transport_views[cell->local_id];
    const auto& outflow_view = cell_outflow_views[cell->local_id];
    const auto& fe_intgrl_values = unit_cell_matrices[cell->local_id];
    const size_t num_nodes = transport_view.GetNumNodes();
    const auto& IntV_shapeI = fe_intgrl_values.intV_shapeI;
    const auto& IntS_shapeI = fe_intgrl_values.intS_shapeI;
    const auto& xs_data = block_id_to_xs_map.at(cell->block_id);

    // Inflow: This is essentially an integration over all faces, all angles, and all groups. For
    // non-reflective boundaries, only the cosines that are negative are added to the inflow
    // integral. For reflective boundaries, it is expected that, upon convergence, inflow = outflow
    // (within numerical tolerances set by the user).
    for (int f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];

      if (not face.has_neighbor) // Boundary face
      {
        const auto& bndry = sweep_boundaries.at(face.neighbor_id);

        if (bndry->IsReflecting())
        {
          for (unsigned int g = 0; g < num_groups; ++g)
            local_in_flow += outflow_view.Get(f, g);
        }
        else
        {
          for (const auto& groupset : groupsets)
          {
            for (size_t n = 0; n < groupset.quadrature->GetNumAngles(); ++n)
            {
              const auto& omega = groupset.quadrature->GetOmega(n);
              const double wt = groupset.quadrature->GetWeight(n);
              const double mu = omega.Dot(face.normal);

              if (mu < 0.0)
              {
                for (int fi = 0; fi < face.vertex_ids.size(); ++fi)
                {
                  const int i = cell_mapping.MapFaceNode(f, fi);
                  const auto& IntFi_shapeI = IntS_shapeI[f](i);

                  for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
                  {
                    const double psi =
                      *bndry->PsiIncoming(cell->local_id, f, fi, n, groupset.id, g);
                    local_in_flow -= mu * wt * psi * IntFi_shapeI;
                  } // for group
                } // for fi
              } // if mu < 0
            } // for n
          } // for groupset
        } // if reflecting boundary
      } // if boundary
    } // for f

    // Outflow: Only physical boundary leakage contributes to the global particle balance.
    for (size_t f = 0; f < cell->faces.size(); ++f)
      if (not cell->faces[f].has_neighbor)
        for (unsigned int g = 0; g < num_groups; ++g)
          local_out_flow += outflow_view.Get(f, g);

    // Absorption and sources
    const auto& xs = transport_view.GetXS();
    const auto& sigma_a = xs.GetSigmaAbsorption();
    const auto& inv_vel = xs.GetInverseVelocity();
    const auto& stopping_power = xs_data->GetStoppingPower();
    const auto charged_ranges = stopping_power.empty()
                                  ? std::vector<std::pair<unsigned int, unsigned int>>{}
                                  : FindChargedGroupRangesForBalance(stopping_power);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (unsigned int g = 0; g < num_groups; ++g)
      {
        auto imap = transport_view.MapDOF(i, 0, g);
        double phi_0g = phi_new_local[imap];
        double q_0g = mat_src[imap];
        if (do_problem.HasUncollidedFlux())
          q_0g -= first_collision_source[imap];

        const double abs_contrib = sigma_a[g] * phi_0g * IntV_shapeI(i);
        local_absorption += abs_contrib;
        local_production += q_0g * IntV_shapeI(i);
      } // for g
    } // for i

    if (options.csda_enabled)
    {
      if (not stopping_power.empty())
      {
        const auto delta_e = xs_data->GetDeltaE();
        const auto& energy_bounds = xs_data->GetEnergyBounds();
        double cell_volume = 0.0;
        for (size_t i = 0; i < num_nodes; ++i)
          cell_volume += IntV_shapeI(i);

        if (cell_volume > 0.0)
        {
          const size_t cell_g_offset = cell->local_id * static_cast<size_t>(num_groups);
          for (size_t r = 0; r < charged_ranges.size(); ++r)
          {
            const auto [g_begin, g_end] = charged_ranges[r];
            const int sign = (r == 0) ? +1 : -1;
            for (unsigned int g = g_begin; g < g_end; ++g)
            {
              double phi_integral = 0.0;
              for (size_t i = 0; i < num_nodes; ++i)
              {
                const auto imap = transport_view.MapDOF(i, 0, g);
                phi_integral += IntV_shapeI(i) * phi_new_local[imap];
              }
              const double phi_avg = phi_integral / cell_volume;
              const double Sg = stopping_power[g];

              local_csda_energy_deposition += Sg * phi_integral;

              if (g + 1 == g_end)
              {
                const double terminal_charge_current =
                  Sg * (phi_avg / delta_e[g] - phi_e_new_local[cell_g_offset + g]);
                local_csda_energy_deposition +=
                  energy_bounds[g + 1] * terminal_charge_current * cell_volume;
                local_csda_charge_deposition += sign * terminal_charge_current * cell_volume;
              }
            }
          }
        }
      }
    }

    if (time_dependent)
    {
      for (const auto& groupset : groupsets)
      {
        const auto& quad = groupset.quadrature;
        const auto num_angles = quad->GetNumAngles();
        const auto first_grp = groupset.first_group;
        const auto num_gs_groups = groupset.GetNumGroups();
        const size_t groupset_angle_group_stride =
          groupset.psi_uk_man_.GetNumberOfUnknowns() * num_gs_groups;
        const size_t groupset_group_stride = num_gs_groups;
        const size_t base = discretization.MapDOFLocal(*cell, 0, groupset.psi_uk_man_, 0, 0);

        for (size_t i = 0; i < num_nodes; ++i)
        {
          for (size_t gsg = 0; gsg < num_gs_groups; ++gsg)
          {
            double phi_old = 0.0;
            double phi_new = 0.0;
            for (size_t n = 0; n < num_angles; ++n)
            {
              const size_t imap =
                base + i * groupset_angle_group_stride + n * groupset_group_stride;
              const double wt = quad->GetWeight(n);
              phi_old += wt * psi_old_local[groupset.id][imap + gsg];
              phi_new += wt * psi_new_local[groupset.id][imap + gsg];
            }

            const auto g = first_grp + gsg;
            const double val = inv_vel[g] * IntV_shapeI(i);
            local_initial += val * phi_old;
            local_final += val * phi_new;
          }
        }
      }
    }
  } // for cell

  // Compute local balance
  double local_balance = local_production + local_in_flow - local_absorption - local_out_flow;
  std::vector<double> local_balance_table = {
    local_absorption, local_production, local_in_flow, local_out_flow, local_balance};
  if (time_dependent)
  {
    local_balance_table.push_back(local_initial);
    local_balance_table.push_back(local_final);
  }
  auto table_size = static_cast<int>(local_balance_table.size());

  // Compute global balance
  std::vector<double> global_balance_table(table_size, 0.0);
  mpi_comm.all_reduce(
    local_balance_table.data(), table_size, global_balance_table.data(), mpi::op::sum<double>());
  double global_absorption = global_balance_table.at(0);
  double global_production = global_balance_table.at(1);
  double global_in_flow = global_balance_table.at(2);
  double global_out_flow = global_balance_table.at(3);
  double global_balance = global_balance_table.at(4);
  if (do_problem.HasUncollidedFlux())
  {
    global_production += do_problem.GetUncollidedSourceRate();
    global_out_flow += do_problem.GetUncollidedOutflowRate();
    global_balance = global_production + global_in_flow - (global_absorption + global_out_flow);
  }
  double global_initial = 0.0;
  double global_final = 0.0;
  double global_csda_charge_deposition = 0.0;
  double global_csda_energy_deposition = 0.0;
  mpi_comm.all_reduce(
    &local_csda_charge_deposition, 1, &global_csda_charge_deposition, mpi::op::sum<double>());
  mpi_comm.all_reduce(
    &local_csda_energy_deposition, 1, &global_csda_energy_deposition, mpi::op::sum<double>());
  double global_csda_particle_balance =
    global_production + global_in_flow -
    (global_absorption + global_out_flow + global_csda_charge_deposition);
  if (scaling_factor != 1.0)
  {
    global_production *= scaling_factor;
    global_balance = global_production + global_in_flow - (global_absorption + global_out_flow);
    global_csda_particle_balance =
      global_production + global_in_flow -
      (global_absorption + global_out_flow + global_csda_charge_deposition);
  }

  if (time_dependent)
  {
    global_initial = global_balance_table.at(5);
    global_final = global_balance_table.at(6);
  }

  do_problem.SetQMomentsFrom(saved_q_moments_local);

  if (time_dependent)
  {
    // NOTE: The time-dependent balance calculation uses end-of-step rates and a simple
    // initial/final calculation. For full consistency with the time discretization, it
    // should use theta-centered rates (sources/absorption), time-averaged boundary
    // inflow/outflow, and a correct handling of pulsed sources/boundaries when
    // start/end times fall inside the timestep.
    const double predicted_inventory_change = dt * global_balance;
    const double actual_inventory_change = global_final - global_initial;
    const double inventory_residual = actual_inventory_change - predicted_inventory_change;
    const BalanceTable table{global_absorption,
                             global_production,
                             global_in_flow,
                             global_out_flow,
                             global_balance,
                             global_csda_charge_deposition,
                             global_csda_particle_balance,
                             global_csda_energy_deposition,
                             global_initial,
                             global_final,
                             predicted_inventory_change,
                             actual_inventory_change,
                             inventory_residual};
    log.Log() << "\nTimestep balance (dt = " << std::setprecision(6) << std::fixed << dt << "):\n"
              << std::setprecision(6) << std::scientific
              << " Production rate             = " << table.production_rate << "\n"
              << " In-flow rate                = " << table.inflow_rate << "\n"
              << " Absorption rate             = " << table.absorption_rate << "\n"
              << " Out-flow rate               = " << table.outflow_rate << "\n"
              << " Balance                     = " << table.balance << "\n"
              << " CSDA charge deposition      = "
              << table.csda_charge_deposition_rate.value_or(0.0) << "\n"
              << " CSDA particle balance       = "
              << table.csda_particle_balance.value_or(0.0) << "\n"
              << " CSDA energy deposition      = "
              << table.csda_energy_deposition_rate.value_or(0.0) << "\n"
              << " Initial inventory           = " << table.initial_inventory.value() << "\n"
              << " Final inventory             = " << table.final_inventory.value() << "\n"
              << " Predicted inventory change  = " << table.predicted_inventory_change.value()
              << "\n"
              << " Actual inventory change     = " << table.actual_inventory_change.value() << "\n"
              << " Inventory residual          = " << table.inventory_residual.value() << "\n\n";

    opensn::mpi_comm.barrier();
    return table;
  }

  opensn::mpi_comm.barrier();
  return {global_absorption,
          global_production,
          global_in_flow,
          global_out_flow,
          global_balance,
          global_csda_charge_deposition,
          global_csda_particle_balance,
          global_csda_energy_deposition,
          std::nullopt,
          std::nullopt,
          std::nullopt,
          std::nullopt,
          std::nullopt};
}

void
ComputeBalance(DiscreteOrdinatesProblem& do_problem, double scaling_factor)
{
  CaliperPhaseScope cali_post_phase("Post", CaliperPostPhaseDepth());
  CALI_CXX_MARK_SCOPE("Balance");

  const auto table = ComputeBalanceTable(do_problem, scaling_factor);
  const auto time_dependent = do_problem.IsTimeDependent();
  const auto dt = time_dependent ? do_problem.GetTimeStep() : 0.0;

  if (time_dependent)
  {
    log.Log() << "\nTimestep balance (dt = " << std::setprecision(6) << std::fixed << dt << "):\n"
              << std::setprecision(6) << std::scientific
              << " Production rate             = " << table.production_rate << "\n"
              << " In-flow rate                = " << table.inflow_rate << "\n"
              << " Absorption rate             = " << table.absorption_rate << "\n"
              << " Out-flow rate               = " << table.outflow_rate << "\n"
              << " Balance                     = " << table.balance << "\n"
              << " CSDA charge deposition      = "
              << table.csda_charge_deposition_rate.value_or(0.0) << "\n"
              << " CSDA particle balance       = "
              << table.csda_particle_balance.value_or(0.0) << "\n"
              << " CSDA energy deposition      = "
              << table.csda_energy_deposition_rate.value_or(0.0) << "\n"
              << " Initial inventory           = " << table.initial_inventory.value_or(0.0) << "\n"
              << " Final inventory             = " << table.final_inventory.value_or(0.0) << "\n"
              << " Predicted inventory change  = " << table.predicted_inventory_change.value_or(0.0)
              << "\n"
              << " Actual inventory change     = " << table.actual_inventory_change.value_or(0.0)
              << "\n"
              << " Inventory residual          = " << table.inventory_residual.value_or(0.0)
              << "\n\n";
  }
  else
  {
    log.Log() << "\nBalance table:\n"
              << std::setprecision(6) << std::scientific
              << " Absorption rate             = " << table.absorption_rate << "\n"
              << " Production rate             = " << table.production_rate << "\n"
              << " In-flow rate                = " << table.inflow_rate << "\n"
              << " Out-flow rate               = " << table.outflow_rate << "\n"
              << " Balance                     = " << table.balance << "\n"
              << " CSDA charge deposition      = "
              << table.csda_charge_deposition_rate.value_or(0.0) << "\n"
              << " CSDA particle balance       = "
              << table.csda_particle_balance.value_or(0.0) << "\n"
              << " CSDA energy deposition      = "
              << table.csda_energy_deposition_rate.value_or(0.0) << "\n\n";
  }
}

std::vector<double>
ComputeLeakage(DiscreteOrdinatesProblem& do_problem,
               const unsigned int groupset_id,
               const uint64_t boundary_id)
{

  OpenSnInvalidArgumentIf(groupset_id >= do_problem.GetNumGroupsets(), "Invalid groupset id.");

  const auto& grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& groupset = do_problem.GetGroupset(groupset_id);
  const auto& cell_transport_views = do_problem.GetCellTransportViews();
  const auto& cell_outflow_views = do_problem.GetCellOutflowViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& sweep_boundaries = do_problem.GetSweepBoundaries();
  const auto& quad = groupset.quadrature;
  const auto num_gs_groups = groupset.GetNumGroups();
  const auto gsi = groupset.first_group;

  std::vector<double> local_leakage(num_gs_groups, 0.0);
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& transport_view = cell_transport_views[cell->local_id];
    const auto& outflow_view = cell_outflow_views[cell->local_id];
    const auto& cell_mapping = sdm.GetCellMapping(*cell);

    for (unsigned int f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];

      if (not face.has_neighbor && face.neighbor_id == boundary_id)
      {
        auto& bndry = *sweep_boundaries.at(face.neighbor_id);

        for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
        {
          auto g = gsi + gsg;
          local_leakage[gsg] += outflow_view.Get(f, g);
          local_leakage[gsg] += GetInflow(
            *cell, cell_mapping, unit_cell_matrices, face, quad, bndry, f, groupset.id, gsg);
        }
      }
    }
  }

  // Reduce
  std::vector<double> global_leakage(num_gs_groups, 0.0);
  mpi_comm.all_reduce(local_leakage.data(),
                      static_cast<int>(num_gs_groups),
                      global_leakage.data(),
                      mpi::op::sum<double>());

  return global_leakage;
}

std::map<uint64_t, std::vector<double>>
ComputeLeakage(DiscreteOrdinatesProblem& do_problem, const std::vector<uint64_t>& boundary_ids)
{

  const auto& grid = do_problem.GetGrid();
  const auto uniq = grid->GetUniqueBoundaryIDs();
  for (const auto& bid : boundary_ids)
    OpenSnInvalidArgumentIf(std::find(uniq.begin(), uniq.end(), bid) == uniq.end(),
                            "Boundary id " + std::to_string(bid) + " not found on grid.");

  const auto num_groups = do_problem.GetNumGroups();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();
  const auto& cell_transport_views = do_problem.GetCellTransportViews();
  const auto& cell_outflow_views = do_problem.GetCellOutflowViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& sweep_boundaries = do_problem.GetSweepBoundaries();

  std::map<uint64_t, std::vector<double>> local_leakage;
  for (const auto& bid : boundary_ids)
    local_leakage[bid].assign(num_groups, 0.0);

  const std::unordered_set<uint64_t> requested_ids(boundary_ids.begin(), boundary_ids.end());

  for (const auto& groupset : groupsets)
  {
    const auto& quad = groupset.quadrature;

    for (const auto& cell : grid->GetLocalCells())
    {
      const auto& transport_view = cell_transport_views[cell->local_id];
      const auto& outflow_view = cell_outflow_views[cell->local_id];
      const auto& cell_mapping = sdm.GetCellMapping(*cell);

      for (unsigned int f = 0; f < cell->faces.size(); ++f)
      {
        const auto& face = cell->faces[f];

        if (not face.has_neighbor && requested_ids.count(face.neighbor_id))
        {
          auto& bndry = *sweep_boundaries.at(face.neighbor_id);

          for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
          {
            auto group_num = groupset.first_group + g;
            local_leakage[face.neighbor_id][group_num] += outflow_view.Get(f, group_num);
            local_leakage[face.neighbor_id][group_num] += GetInflow(
              *cell, cell_mapping, unit_cell_matrices, face, quad, bndry, f, groupset.id, g);
          }
        }
      }
    }
  }

  // Serialize in user-requested order
  std::vector<double> local_data;
  local_data.reserve(boundary_ids.size() * num_groups);
  for (const auto& bid : boundary_ids)
  {
    const auto& vec = local_leakage.at(bid);
    local_data.insert(local_data.end(), vec.begin(), vec.end());
  }

  // Reduce
  auto count = static_cast<const int>(local_data.size());
  std::vector<double> global_data(count, 0.0);
  mpi_comm.all_reduce(local_data.data(), count, global_data.data(), mpi::op::sum<double>());

  // Unpack in user-requested order
  std::map<uint64_t, std::vector<double>> global_leakage;
  for (const auto& bid : boundary_ids)
    global_leakage[bid] = std::vector<double>();

  for (size_t b = 0; b < boundary_ids.size(); ++b)
  {
    auto& vec = global_leakage[boundary_ids[b]];
    vec.reserve(num_groups);
    for (unsigned int g = 0; g < num_groups; ++g)
      vec.push_back(global_data[b * num_groups + g]);
  }

  return global_leakage;
}

} // namespace opensn
