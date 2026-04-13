// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/utils/error.h"
#include <iomanip>
#include <memory>
#include <sstream>

namespace opensn
{

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
    const auto num_angles = groupset.quadrature->omegas.size();

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

  for (const auto& cell : grid_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto imapA = sdm.MapDOFLocal(cell, i, uk_man, angle, group_in_groupset);
      const auto imapB = sdm.MapDOFLocal(cell, i);
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

} // namespace opensn
