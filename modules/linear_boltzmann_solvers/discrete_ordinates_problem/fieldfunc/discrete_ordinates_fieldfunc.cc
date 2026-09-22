// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/csda_utils.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/utils/error.h"
#include <iomanip>
#include <map>
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
  const bool is_cepxs_energy_deposition = xs_name == "cepxs_energy_deposition";
  const bool is_csda_charge_deposition = xs_name == "csda_charge_deposition";
  const bool is_csda_charge_deposition_term = xs_name == "csda_charge_deposition_term";
  const bool is_csda_charge_deposition_term_cellavg =
    xs_name == "csda_charge_deposition_term_cellavg";
  if (not is_csda_energy_deposition and not is_cepxs_energy_deposition and
      not is_csda_charge_deposition and not is_csda_charge_deposition_term and
      not is_csda_charge_deposition_term_cellavg)
    return std::nullopt;

  OpenSnInvalidArgumentIf(not options_.csda_enabled,
                          GetName() + ": Field function \"" + xs_name +
                            "\" requires `options.csda_enabled=true`.");

  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;
  const auto& unit_cell_matrices = GetUnitCellMatrices();
  std::vector<double> data_vector_local(local_node_count_, 0.0);

  // The charge sign depends on the problem-level block a group belongs to, not on the
  // ranges of an individual material.
  const auto problem_charged_ranges =
    FindCSDAProblemChargedGroupRanges(block_id_to_xs_map_, num_groups_);

  const auto energy = MultiGroupXS::ResolveEnergyGroupStructure(block_id_to_xs_map_, num_groups_);

  // Per-material data, built once on the first local cell of each block.
  struct MaterialData
  {
    const MultiGroupXS* xs = nullptr;
    std::vector<double> delta_e;
    std::vector<std::pair<unsigned int, unsigned int>> charged_ranges;
    std::vector<double> conservative_energy_coeffs;
    const std::vector<double>* raw_coeffs = nullptr;
  };
  std::map<unsigned int, MaterialData> material_data;

  const auto GetMaterialData = [&](const unsigned int block_id) -> const MaterialData&
  {
    if (const auto it = material_data.find(block_id); it != material_data.end())
      return it->second;

    auto& data = material_data[block_id];
    const auto& xs = *block_id_to_xs_map_.at(block_id);
    data.xs = &xs;
    const auto& stopping_power = xs.GetStoppingPower();
    if (not stopping_power.empty())
    {
      data.delta_e = energy.widths;
      OpenSnLogicalErrorIf(stopping_power.size() != num_groups_,
                           GetName() + ": Field function \"" + xs_name +
                             "\" requires stopping power with one value per group.");
      OpenSnLogicalErrorIf(data.delta_e.size() != num_groups_,
                           GetName() + ": Field function \"" + xs_name +
                             "\" requires energy-bin widths with one value per group.");
      data.charged_ranges = xs.GetStoppingPowerGroupRanges();
    }

    if (is_csda_energy_deposition)
    {
      data.conservative_energy_coeffs = xs.ComputeCollisionEnergyLossCoefficients(energy);
      data.raw_coeffs = &data.conservative_energy_coeffs;
    }
    else if (is_cepxs_energy_deposition)
      data.raw_coeffs = xs.GetEnergyDeposition().empty() ? nullptr : &xs.GetEnergyDeposition();
    else if (is_csda_charge_deposition)
      data.raw_coeffs =
        xs.HasCustomXS("charge_deposition") ? &xs.GetCustomXS("charge_deposition") : nullptr;

    // A material with no charged groups and no charge-deposition data contributes zero.
    OpenSnLogicalErrorIf(is_csda_charge_deposition and data.raw_coeffs == nullptr and
                           not data.charged_ranges.empty(),
                         GetName() + ": Field function \"" + xs_name +
                           "\" could not find custom XS \"charge_deposition\".");
    if (data.raw_coeffs != nullptr)
      OpenSnLogicalErrorIf(data.raw_coeffs->size() != num_groups_,
                           GetName() + ": Field function \"" + xs_name +
                             "\" has incompatible raw deposition XS group size.");
    return data;
  };

  for (const auto& cell : grid_->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();
    const auto& data = GetMaterialData(cell->block_id);
    const auto& stopping_power = data.xs->GetStoppingPower();

    const size_t cell_g_offset = cell->local_id * static_cast<size_t>(num_groups_);

    if (is_csda_charge_deposition_term)
    {
      // Preserve the normal PWLD nodal representation for this public field. The energy-slope
      // moment is cellwise, while the scalar flux contribution is evaluated at each node.
      for (size_t i = 0; i < num_nodes; ++i)
      {
        double nodal_value = 0.0;
        for (const auto& charged_range : data.charged_ranges)
        {
          const unsigned int g = charged_range.second - 1;
          const auto imap = sdm.MapDOFLocal(*cell, i, phi_uk_man, 0, g);
          const double slowing_down_density =
            stopping_power[g] *
            (phi_new_local_[imap] / data.delta_e[g] - phi_e_new_local_[cell_g_offset + g]);
          nodal_value += CSDAChargeSign(problem_charged_ranges, g) * slowing_down_density;
        }
        data_vector_local[sdm.MapDOFLocal(*cell, i)] = nodal_value;
      }
      continue;
    }

    const auto phi_cell_avg = ComputeCellAveragePhi0g(
      sdm, *cell, phi_uk_man, unit_cell_matrices, phi_new_local_, num_groups_);

    // The CSDA terminal correction below is fundamentally a cell-averaged quantity: it is
    // built entirely from cell-averaged flux and psi_E moments, so it is constant across a cell.
    // The reaction-rate term can vary within a cell when evaluated from the nodal flux. Summing a
    // nodal quantity against the cell-flat CSDA correction can leave a spurious intra-cell
    // variation when the two terms nearly cancel. Evaluate both at cell resolution so their sum
    // represents one consistent discrete energy-deposition quantity.
    double cell_avg_raw_value = 0.0;
    if (data.raw_coeffs != nullptr)
      for (unsigned int g = 0; g < num_groups_; ++g)
        cell_avg_raw_value += (*data.raw_coeffs)[g] * phi_cell_avg[g];

    double csda_value = 0.0;
    for (const auto& [g_begin, g_end] : data.charged_ranges)
    {
      for (unsigned int g = g_begin; g < g_end; ++g)
      {
        const bool is_terminal_group = (g + 1 == g_end);
        // Slowing-down density at the low-energy edge of group g, from the within-group
        // ansatz: psi(E_low) = Phi_g / DeltaE_g - psiE_g.
        const double slowing_down_density =
          stopping_power[g] *
          (phi_cell_avg[g] / data.delta_e[g] - phi_e_new_local_[cell_g_offset + g]);

        if (is_csda_energy_deposition or is_cepxs_energy_deposition)
        {
          // Conservative response: weight each edge by the decrease in midpoint energy.
          const double next_energy = is_terminal_group ? 0.0 : energy.centers[g + 1];
          const double edge_energy_loss = energy.centers[g] - next_energy;
          csda_value += edge_energy_loss * slowing_down_density;
        }
        else if (is_terminal_group)
        {
          // Particles leaving the terminal group are deposited, carrying the charge of
          // their problem-level block.
          csda_value += CSDAChargeSign(problem_charged_ranges, g) * slowing_down_density;
        }
      }
    }

    const double cell_value = cell_avg_raw_value + csda_value;

    for (size_t i = 0; i < num_nodes; ++i)
      data_vector_local[sdm.MapDOFLocal(*cell, i)] = cell_value;
  }

  return data_vector_local;
}

} // namespace opensn
