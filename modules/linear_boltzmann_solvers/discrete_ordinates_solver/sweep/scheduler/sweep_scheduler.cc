// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds_adams_adams_hawkins.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/reflecting_boundary.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <sstream>
#include <algorithm>

namespace opensn
{

SweepScheduler::SweepScheduler(SchedulingAlgorithm scheduler_type,
                               AngleAggregation& angle_agg,
                               SweepChunk& sweep_chunk)
  : scheduler_type_(scheduler_type), angle_agg_(angle_agg), sweep_chunk_(sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::SweepScheduler");

  angle_agg_.InitializeReflectingBCs();

  if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    InitializeAlgoDOG();

  // Initialize delayed upstream data
  for (auto& angsetgrp : angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      angset->InitializeDelayedUpstreamData();

  // Get local max num messages accross anglesets
  int local_max_num_messages = 0;
  for (auto& angsetgrp : angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      local_max_num_messages = std::max(angset->GetMaxBufferMessages(), local_max_num_messages);

  // Reconcile all local maximums
  int global_max_num_messages = 0;
  mpi_comm.all_reduce(local_max_num_messages, global_max_num_messages, mpi::op::max<int>());

  // Propogate items back to sweep buffers
  for (auto& angsetgrp : angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      angset->SetMaxBufferMessages(global_max_num_messages);
}

SweepChunk&
SweepScheduler::GetSweepChunk()
{
  return sweep_chunk_;
}

void
SweepScheduler::InitializeAlgoDOG()
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::InitializeAlgoDOG");

  // Load all anglesets in preperation for sorting
  // Loop over angleset groups
  for (size_t q = 0; q < angle_agg_.angle_set_groups.size(); q++)
  {
    AngleSetGroup& angleset_group = angle_agg_.angle_set_groups[q];

    // Loop over anglesets in group
    size_t num_anglesets = angleset_group.AngleSets().size();
    for (size_t as = 0; as < num_anglesets; as++)
    {
      auto angleset = angleset_group.AngleSets()[as];
      const auto& spds = dynamic_cast<const SPDS_AdamsAdamsHawkins&>(angleset->GetSPDS());

      const std::vector<STDG>& leveled_graph = spds.GetGlobalSweepPlanes();

      // Find location depth
      int loc_depth = -1;
      for (size_t level = 0; level < leveled_graph.size(); level++)
      {
        for (size_t index = 0; index < leveled_graph[level].item_id.size(); index++)
        {
          if (leveled_graph[level].item_id[index] == opensn::mpi_comm.rank())
          {
            loc_depth = static_cast<int>(leveled_graph.size() - level);
            break;
          }
        } // for locations in plane
      }   // for sweep planes

      // Set up rule values
      if (loc_depth >= 0)
      {
        RuleValues new_rule_vals(angleset);
        new_rule_vals.depth_of_graph = loc_depth;
        new_rule_vals.set_index = as + q * num_anglesets;

        const auto& omega = spds.Omega();
        new_rule_vals.sign_of_omegax = (omega.x >= 0) ? 2 : 1;
        new_rule_vals.sign_of_omegay = (omega.y >= 0) ? 2 : 1;
        new_rule_vals.sign_of_omegaz = (omega.z >= 0) ? 2 : 1;

        rule_values_.push_back(new_rule_vals);
      }
      else
      {
        log.LogAllError() << "Location depth not found in Depth-Of-Graph algorithm.";
        Exit(EXIT_FAILURE);
      }

    } // for anglesets
  }   // for quadrants/anglesetgroups

  // Init sort functions
  struct
  {
    bool operator()(const RuleValues& a, const RuleValues& b)
    {
      return a.depth_of_graph > b.depth_of_graph;
    }
  } compare_D;

  struct
  {
    bool operator()(const RuleValues& a, const RuleValues& b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and (a.sign_of_omegax > b.sign_of_omegax);
    }
  } compare_omega_x;

  struct
  {
    bool operator()(const RuleValues& a, const RuleValues& b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and (a.sign_of_omegax == b.sign_of_omegax) and
             (a.sign_of_omegay > b.sign_of_omegay);
    }
  } compare_omega_y;

  struct
  {
    bool operator()(const RuleValues& a, const RuleValues& b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and (a.sign_of_omegax == b.sign_of_omegax) and
             (a.sign_of_omegay == b.sign_of_omegay) and (a.sign_of_omegaz > b.sign_of_omegaz);
    }
  } compare_omega_z;

  // Sort
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_D);
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_omega_x);
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_omega_y);
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_omega_z);
}

void
SweepScheduler::ScheduleAlgoDOG(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoDOG");

  // Loop till done
  bool finished = false;
  while (not finished)
  {
    finished = true;
    for (auto& rule_value : rule_values_)
    {
      auto angleset = rule_value.angle_set;

      // Query angleset status
      // Angleset status will be one of the following:
      //  - RECEIVING.
      //      Meaning it is either waiting for messages or actively receiving it
      //  - READY_TO_EXECUTE.
      //      Meaning it has received all upstream data and can be executed
      //  - FINISHED.
      //      Meaning the angleset has executed its sweep chunk
      AngleSetStatus status =
        angleset->AngleSetAdvance(sweep_chunk, AngleSetStatus::NO_EXEC_IF_READY);

      // Execute if ready and allowed
      // If this angleset is the one scheduled to run
      // and it is ready then it will be given permission
      if (status == AngleSetStatus::READY_TO_EXECUTE)
      {
        std::stringstream message_i;
        message_i << "Angleset " << angleset->GetID() << " executed on location "
                  << opensn::mpi_comm.rank();

        status = angleset->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);

        std::stringstream message_f;
        message_f << "Angleset " << angleset->GetID() << " finished on location "
                  << opensn::mpi_comm.rank();
      }

      if (status != AngleSetStatus::FINISHED)
        finished = false;
    } // for each angleset rule
  }   // while not finished

  // Receive delayed data
  opensn::mpi_comm.barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.AngleSets())
      {
        if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
          received_delayed_data = false;

        if (not angle_set->ReceiveDelayedData())
          received_delayed_data = false;
      }
  }

  // Reset all
  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.AngleSets())
      angle_set->ResetSweepBuffers();

  for (auto& [bid, bndry] : angle_agg_.GetSimBoundaries())
  {
    if (bndry->Type() == BoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<ReflectingBoundary>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }
}

void
SweepScheduler::ScheduleAlgoFIFO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoFIFO");

  // Loop over AngleSetGroups
  AngleSetStatus completion_status = AngleSetStatus::NOT_FINISHED;
  while (completion_status == AngleSetStatus::NOT_FINISHED)
  {
    completion_status = AngleSetStatus::FINISHED;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.AngleSets())
      {
        const auto angle_set_status =
          angle_set->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
        if (angle_set_status == AngleSetStatus::NOT_FINISHED)
          completion_status = AngleSetStatus::NOT_FINISHED;
      } // for angleset
  }     // while not finished

  // Receive delayed data
  opensn::mpi_comm.barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.AngleSets())
      {
        if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
          received_delayed_data = false;

        if (not angle_set->ReceiveDelayedData())
          received_delayed_data = false;
      }
  }

  // Reset all
  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.AngleSets())
      angle_set->ResetSweepBuffers();

  for (auto& [bid, bndry] : angle_agg_.GetSimBoundaries())
  {
    if (bndry->Type() == BoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<ReflectingBoundary>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }
}

void
SweepScheduler::Sweep()
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::Sweep");

  if (scheduler_type_ == SchedulingAlgorithm::FIRST_IN_FIRST_OUT)
    ScheduleAlgoFIFO(sweep_chunk_);
  else if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    ScheduleAlgoDOG(sweep_chunk_);
}

void
SweepScheduler::SetDestinationPhi(std::vector<double>& destination_phi)
{
  sweep_chunk_.SetDestinationPhi(destination_phi);
}

void
SweepScheduler::ZeroDestinationPhi()
{
  sweep_chunk_.ZeroDestinationPhi();
}

std::vector<double>&
SweepScheduler::GetDestinationPhi()
{
  return sweep_chunk_.GetDestinationPhi();
}

void
SweepScheduler::SetDestinationPsi(std::vector<double>& destination_psi)
{
  sweep_chunk_.SetDestinationPsi(destination_psi);
}

void
SweepScheduler::ZeroDestinationPsi()
{
  sweep_chunk_.ZeroDestinationPsi();
}

std::vector<double>&
SweepScheduler::GetDestinationPsi()
{
  return sweep_chunk_.GetDestinationPsi();
}

void
SweepScheduler::ZeroIncomingDelayedPsi()
{
  angle_agg_.ZeroIncomingDelayedPsi();
}

void
SweepScheduler::ZeroOutgoingDelayedPsi()
{
  angle_agg_.ZeroOutgoingDelayedPsi();
}

void
SweepScheduler::ZeroOutputFluxDataStructures()
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ZeroOutputFluxDataStructures");

  ZeroDestinationPsi();
  ZeroDestinationPhi();
  ZeroOutgoingDelayedPsi();
}

void
SweepScheduler::SetBoundarySourceActiveFlag(bool flag_value)
{
  sweep_chunk_.SetBoundarySourceActiveFlag(flag_value);
}

} // namespace opensn
