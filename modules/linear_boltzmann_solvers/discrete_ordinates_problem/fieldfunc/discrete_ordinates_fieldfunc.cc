// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/csda_utils.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/utils/error.h"
#include <iomanip>
#include <memory>
#include <optional>
#include <sstream>

namespace opensn
{
namespace
{

std::vector<double>
ComputeCellAveragePhi0g(const SpatialDiscretization& sdm,
                        const Cell& cell,
                        const UnknownManager& phi_uk_man,
                        const std::vector<UnitCellMatrices>& unit_cell_matrices,
                        const std::vector<double>& phi_new_local,
                        const unsigned int num_groups)
{
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.GetNumNodes();
  const auto& intV_shapeI = unit_cell_matrices[cell.local_id].intV_shapeI;

  double cell_volume = 0.0;
  for (size_t i = 0; i < num_nodes; ++i)
    cell_volume += intV_shapeI(i);

  std::vector<double> phi_avg(num_groups, 0.0);
  if (cell_volume <= 0.0)
    return phi_avg;

  for (unsigned int g = 0; g < num_groups; ++g)
  {
    double phi_integral = 0.0;
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto imap = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, g);
      phi_integral += intV_shapeI(i) * phi_new_local[imap];
    }
    phi_avg[g] = phi_integral / cell_volume;
  }

  return phi_avg;
}
} // namespace

std::vector<std::shared_ptr<FieldFunctionGridBased>>
DiscreteOrdinatesProblem::CreateAngularFluxFieldFunctionList(
  const std::vector<unsigned int>& groups, const std::vector<size_t>& angles)
{
  OpenSnLogicalErrorIf(discretization_ == nullptr || grid_ == nullptr || groupsets_.empty(),
                       "CreateAngularFluxFieldFunctionList: problem not fully constructed.");

  OpenSnLogicalErrorIf(groups.empty(),
                       "CreateAngularFluxFieldFunctionList: groups cannot be empty.");
  OpenSnLogicalErrorIf(angles.empty(),
                       "CreateAngularFluxFieldFunctionList: angles cannot be empty.");

  std::vector<std::shared_ptr<FieldFunctionGridBased>> result;
  result.reserve(groups.size() * angles.size());

  for (const auto g : groups)
  {
    OpenSnLogicalErrorIf(g >= num_groups_,
                         "CreateAngularFluxFieldFunctionList: group index out of range.");

    const LBSGroupset* gs_ptr = nullptr;
    size_t gs_id = 0;
    for (const auto& groupset : groupsets_)
      if (g >= groupset.first_group && g <= groupset.last_group)
      {
        gs_ptr = &groupset;
        gs_id = groupset.id;
        break;
      }
    OpenSnLogicalErrorIf(gs_ptr == nullptr,
                         "CreateAngularFluxFieldFunctionList: group not found in any groupset.");

    const auto& groupset = *gs_ptr;
    const auto num_angles = groupset.quadrature->GetNumAngles();

    for (const auto a : angles)
    {
      OpenSnLogicalErrorIf(a >= num_angles,
                           "CreateAngularFluxFieldFunctionList: angle index out of range for "
                           "groupset " +
                             std::to_string(gs_id) + ".");

      auto ff_ptr = CreateEmptyFieldFunction(MakeAngularFieldFunctionName(gs_id, g, a));
      UpdateAngularFluxFieldFunction(*ff_ptr, gs_id, g, a);
      const std::weak_ptr<LBSProblem> weak_owner = weak_from_this();
      ff_ptr->SetUpdateCallback(
        [weak_owner, gs_id, g, a](FieldFunctionGridBased& ff)
        {
          auto owner = weak_owner.lock();
          OpenSnLogicalErrorIf(not owner,
                               "Cannot update field function after its owning problem has "
                               "been destroyed.");
          auto do_owner = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(owner);
          OpenSnLogicalErrorIf(not do_owner,
                               "Angular flux field function owner is not a "
                               "DiscreteOrdinatesProblem.");
          do_owner->UpdateAngularFluxFieldFunction(ff, gs_id, g, a);
        },
        [weak_owner]() { return not weak_owner.expired(); });
      result.push_back(ff_ptr);
    }
  }

  return result;
}

std::string
DiscreteOrdinatesProblem::MakeAngularFieldFunctionName(const size_t groupset_id,
                                                       const unsigned int group,
                                                       const size_t angle) const
{
  std::ostringstream oss;
  oss << MakeFieldFunctionName("psi_g") << std::setw(3) << std::setfill('0')
      << static_cast<int>(group) << "_a" << std::setw(3) << std::setfill('0')
      << static_cast<int>(angle) << "_gs" << std::setw(2) << std::setfill('0')
      << static_cast<int>(groupset_id);
  return oss.str();
}

std::vector<double>
DiscreteOrdinatesProblem::ComputeAngularFieldFunctionData(const size_t groupset_id,
                                                          const unsigned int group,
                                                          const size_t angle) const
{
  std::vector<double> data_vector_local(local_node_count_, 0.0);

  if (groupset_id >= psi_new_local_.size() || psi_new_local_[groupset_id].empty())
    return data_vector_local;

  const auto& sdm = *discretization_;
  const auto& groupset = groupsets_.at(groupset_id);
  const auto group_start = static_cast<size_t>(groupset.first_group);
  const auto group_in_groupset = group - group_start;
  const auto& uk_man = groupset.psi_uk_man_;
  const auto& psi = psi_new_local_.at(groupset_id);

  for (const auto& cell : grid_->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto imapA = sdm.MapDOFLocal(*cell, i, uk_man, angle, group_in_groupset);
      const auto imapB = sdm.MapDOFLocal(*cell, i);
      data_vector_local[imapB] = psi[imapA];
    }
  }

  return data_vector_local;
}

void
DiscreteOrdinatesProblem::UpdateAngularFluxFieldFunction(FieldFunctionGridBased& ff,
                                                         const size_t groupset_id,
                                                         const unsigned int group,
                                                         const size_t angle)
{
  ff.UpdateFieldVector(ComputeAngularFieldFunctionData(groupset_id, group, angle));
}

std::optional<std::vector<double>>
DiscreteOrdinatesProblem::ComputeDerivedFieldFunctionData(const std::string& xs_name) const
{
  const bool is_energy_deposition_alias = options_.csda_enabled and xs_name == "energy_deposition";
  const bool is_csda_energy_deposition =
    xs_name == "csda_energy_deposition" or is_energy_deposition_alias;
  const bool is_csda_charge_deposition = xs_name == "csda_charge_deposition";
  const bool is_csda_charge_deposition_term = xs_name == "csda_charge_deposition_term";
  const bool is_csda_charge_deposition_term_cellavg =
    xs_name == "csda_charge_deposition_term_cellavg";
  if (not is_csda_energy_deposition and not is_csda_charge_deposition and
      not is_csda_charge_deposition_term and not is_csda_charge_deposition_term_cellavg)
    return std::nullopt;

  OpenSnInvalidArgumentIf(not options_.csda_enabled,
                          GetName() + ": Field function \"" + xs_name +
                            "\" requires `options.csda_enabled=true`.");

  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;
  const auto& unit_cell_matrices = GetUnitCellMatrices();
  std::vector<double> data_vector_local(local_node_count_, 0.0);

  for (const auto& cell : grid_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto& xs = block_id_to_xs_map_.at(cell.block_id);
    const auto& stopping_power = xs->GetStoppingPower();
    const auto delta_e = stopping_power.empty() ? std::vector<double>{} : xs->GetDeltaE();
    const auto& energy_bounds = xs->GetEnergyBounds();

    if (not stopping_power.empty())
    {
      OpenSnLogicalErrorIf(stopping_power.size() != num_groups_,
                           GetName() + ": Field function \"" + xs_name +
                             "\" requires stopping power with one value per group.");
      OpenSnLogicalErrorIf(delta_e.size() != num_groups_,
                           GetName() + ": Field function \"" + xs_name +
                             "\" requires energy-bin widths with one value per group.");
      OpenSnLogicalErrorIf(energy_bounds.size() != num_groups_ + 1,
                           GetName() + ": Field function \"" + xs_name +
                             "\" requires energy-bin boundaries.");
    }

    const std::vector<double>* raw_coeffs = nullptr;
    if (is_csda_energy_deposition)
      raw_coeffs = xs->GetEnergyDeposition().empty() ? nullptr : &xs->GetEnergyDeposition();
    else if (is_csda_charge_deposition)
      raw_coeffs =
        xs->HasCustomXS("charge_deposition") ? &xs->GetCustomXS("charge_deposition") : nullptr;

    OpenSnLogicalErrorIf(is_csda_charge_deposition and raw_coeffs == nullptr,
                         GetName() + ": Field function \"" + xs_name +
                           "\" could not find "
                           "custom XS \"charge_deposition\".");

    if (raw_coeffs != nullptr)
      OpenSnLogicalErrorIf(raw_coeffs->size() != num_groups_,
                           GetName() + ": Field function \"" + xs_name +
                             "\" has incompatible raw deposition XS group size.");

    const auto charged_ranges = stopping_power.empty()
                                  ? std::vector<std::pair<unsigned int, unsigned int>>{}
                                  : FindCSDAChargedGroupRanges(stopping_power);
    if (is_csda_charge_deposition)
      OpenSnLogicalErrorIf(charged_ranges.size() > 2,
                           GetName() + ": Field function \"" + xs_name +
                             "\" supports at most two charged-particle blocks "
                             "(electron then positron).");

    const size_t cell_g_offset = cell.local_id * static_cast<size_t>(num_groups_);
    const auto phi_cell_avg = ComputeCellAveragePhi0g(
      sdm, cell, phi_uk_man, unit_cell_matrices, phi_new_local_, num_groups_);

    // The CSDA terminal-current correction below is fundamentally a cell-averaged quantity: it is
    // built entirely from cell-averaged flux and psi_E moments, so it is constant across a cell.
    // The raw reaction-rate term (from a custom "energy_deposition"/"charge_deposition" XS) is
    // instead evaluated from the true nodal flux and can vary meaningfully within a cell. Summing
    // a nodal quantity against a cell-flat one is harmless when one of the two is negligible, but
    // for CEPXS-CSDA libraries where the raw XS is populated within the charged-particle block
    // (e.g. coupled electron-photon libraries), the two terms can be the same order of magnitude
    // and nearly cancel. Any real intra-cell slope in the nodal term then survives uncanceled and
    // is amplified relative to the much smaller combined signal, producing a spurious sawtooth
    // that has nothing to do with the transport solution itself. Cell-averaging the raw term
    // before combining -- i.e. representing both terms at the same (cell) resolution -- removes
    // that artifact without changing the mean value of either field.
    double cell_avg_raw_value = 0.0;
    if (raw_coeffs != nullptr)
      for (unsigned int g = 0; g < num_groups_; ++g)
        cell_avg_raw_value += raw_coeffs->at(g) * phi_cell_avg[g];

    double csda_value = 0.0;
    for (size_t r = 0; r < charged_ranges.size(); ++r)
    {
      const auto [g_begin, g_end] = charged_ranges[r];
      for (unsigned int g = g_begin; g < g_end; ++g)
      {
        const double phi_g = phi_cell_avg[g];
        const double phi_e_g = phi_e_new_local_[cell_g_offset + g];
        const double Sg = stopping_power[g];
        const double dEg = delta_e[g];
        const bool is_terminal_group = (g + 1 == g_end);
        const double terminal_charge_current = Sg * (phi_g / dEg - phi_e_g);
        const double next_energy =
          is_terminal_group ? 0.0 : CSDAGroupCenterEnergy(energy_bounds, g + 1);
        const double edge_energy_loss = CSDAGroupCenterEnergy(energy_bounds, g) - next_energy;

        if (is_csda_energy_deposition)
        {
          // Conservative response for the discrete CSDA energy-space current.
          csda_value += edge_energy_loss * terminal_charge_current;
        }
        else
        {
          // Derived from the OpenSn within-group ansatz evaluated at the low-energy
          // edge of the terminal charged group:
          //   psi(E_low) = Phi_g / DeltaE_g - psiE_g.
          // Charge deposition is the terminal energy-space current S_g * psi(E_low),
          // positive for electrons and negative for positrons.
          if (is_terminal_group)
            csda_value += ComputeCSDATerminalDepositionRates(terminal_charge_current, r).charge;
        }
      }
    }

    const double cell_value =
      is_csda_charge_deposition_term_cellavg ? csda_value : cell_avg_raw_value + csda_value;

    for (size_t i = 0; i < num_nodes; ++i)
      data_vector_local[sdm.MapDOFLocal(cell, i)] = cell_value;
  }

  return data_vector_local;
}

} // namespace opensn
