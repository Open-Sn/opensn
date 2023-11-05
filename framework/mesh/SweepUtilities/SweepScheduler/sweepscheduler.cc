#include "opensn/framework/mesh/SweepUtilities/SweepScheduler/sweepscheduler.h"
#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "opensn/framework/mesh/SweepUtilities/SPDS/SPDS_AdamsAdamsHawkins.h"
#include "opensn/framework/mesh/SweepUtilities/SweepBoundary/boundary_reflecting.h"
#include "opensn/framework/mpi/chi_mpi.h"
#include <sstream>
#include <algorithm>

namespace chi_mesh::sweep_management
{

SweepScheduler::SweepScheduler(SchedulingAlgorithm in_scheduler_type,
                               AngleAggregation& in_angle_agg,
                               SweepChunk& in_sweep_chunk)
  : scheduler_type_(in_scheduler_type),
    angle_agg_(in_angle_agg),
    sweep_chunk_(in_sweep_chunk),
    sweep_event_tag_(Chi::log.GetRepeatingEventTag("Sweep Timing")),
    sweep_timing_events_tag_(
      {Chi::log.GetRepeatingEventTag("Sweep Chunk Only Timing"), sweep_event_tag_})

{
  angle_agg_.InitializeReflectingBCs();

  if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH) InitializeAlgoDOG();

  //=================================== Initialize delayed upstream data
  for (auto& angsetgrp : in_angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      angset->InitializeDelayedUpstreamData();

  //=================================== Get local max num messages accross
  //                                    anglesets
  int local_max_num_messages = 0;
  for (auto& angsetgrp : in_angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      local_max_num_messages = std::max(angset->GetMaxBufferMessages(), local_max_num_messages);

  //=================================== Reconcile all local maximums
  int global_max_num_messages = 0;
  MPI_Allreduce(
    &local_max_num_messages, &global_max_num_messages, 1, MPI_INT, MPI_MAX, Chi::mpi.comm);

  //=================================== Propogate items back to sweep buffers
  for (auto& angsetgrp : in_angle_agg.angle_set_groups)
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
  //================================================== Load all anglesets
  //                                                   in preperation for
  //                                                   sorting
  //======================================== Loop over angleset groups
  for (size_t q = 0; q < angle_agg_.angle_set_groups.size(); q++)
  {
    TAngleSetGroup& angleset_group = angle_agg_.angle_set_groups[q];

    //================================= Loop over anglesets in group
    size_t num_anglesets = angleset_group.AngleSets().size();
    for (size_t as = 0; as < num_anglesets; as++)
    {
      auto angleset = angleset_group.AngleSets()[as];
      const auto& spds = dynamic_cast<const SPDS_AdamsAdamsHawkins&>(angleset->GetSPDS());

      const TLEVELED_GRAPH& leveled_graph = spds.GetGlobalSweepPlanes();

      //========================== Find location depth
      int loc_depth = -1;
      for (size_t level = 0; level < leveled_graph.size(); level++)
      {
        for (size_t index = 0; index < leveled_graph[level].item_id.size(); index++)
        {
          if (leveled_graph[level].item_id[index] == Chi::mpi.location_id)
          {
            loc_depth = static_cast<int>(leveled_graph.size() - level);
            break;
          }
        } // for locations in plane
      }   // for sweep planes

      //========================== Set up rule values
      if (loc_depth >= 0)
      {
        RULE_VALUES new_rule_vals(angleset);
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
        Chi::log.LogAllError() << "Location depth not found in Depth-Of-Graph algorithm.";
        Chi::Exit(EXIT_FAILURE);
      }

    } // for anglesets
  }   // for quadrants/anglesetgroups

  //================================================== Init sort functions
  struct
  {
    bool operator()(const RULE_VALUES& a, const RULE_VALUES& b)
    {
      return a.depth_of_graph > b.depth_of_graph;
    }
  } compare_D;

  struct
  {
    bool operator()(const RULE_VALUES& a, const RULE_VALUES& b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and (a.sign_of_omegax > b.sign_of_omegax);
    }
  } compare_omega_x;

  struct
  {
    bool operator()(const RULE_VALUES& a, const RULE_VALUES& b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and (a.sign_of_omegax == b.sign_of_omegax) and
             (a.sign_of_omegay > b.sign_of_omegay);
    }
  } compare_omega_y;

  struct
  {
    bool operator()(const RULE_VALUES& a, const RULE_VALUES& b)
    {
      return (a.depth_of_graph == b.depth_of_graph) and (a.sign_of_omegax == b.sign_of_omegax) and
             (a.sign_of_omegay == b.sign_of_omegay) and (a.sign_of_omegaz > b.sign_of_omegaz);
    }
  } compare_omega_z;

  //================================================== Sort
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_D);
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_omega_x);
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_omega_y);
  std::stable_sort(rule_values_.begin(), rule_values_.end(), compare_omega_z);
}

void
SweepScheduler::ScheduleAlgoDOG(SweepChunk& sweep_chunk)
{
  typedef ExecutionPermission ExePerm;
  typedef AngleSetStatus Status;

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::EVENT_BEGIN);

  auto ev_info = std::make_shared<chi::ChiLog::EventInfo>(std::string("Sweep initiated"));

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::SINGLE_OCCURRENCE, ev_info);

  //==================================================== Loop till done
  bool finished = false;
  size_t scheduled_angleset = 0;
  while (!finished)
  {
    finished = true;
    for (auto& rule_value : rule_values_)
    {
      auto angleset = rule_value.angle_set;

      //=============================== Query angleset status
      // Status will here be one of the following:
      //  - RECEIVING.
      //      Meaning it is either waiting for messages or actively receiving it
      //  - READY_TO_EXECUTE.
      //      Meaning it has received all upstream data and can be executed
      //  - FINISHED.
      //      Meaning the angleset has executed its sweep chunk
      Status status =
        angleset->AngleSetAdvance(sweep_chunk, sweep_timing_events_tag_, ExePerm::NO_EXEC_IF_READY);

      //=============================== Execute if ready and allowed
      // If this angleset is the one scheduled to run
      // and it is ready then it will be given permission
      if (status == Status::READY_TO_EXECUTE)
      {
        std::stringstream message_i;
        message_i << "Angleset " << angleset->GetID() << " executed on location "
                  << Chi::mpi.location_id;

        auto ev_info_i = std::make_shared<chi::ChiLog::EventInfo>(message_i.str());

        Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::SINGLE_OCCURRENCE, ev_info_i);

        status = angleset->AngleSetAdvance(sweep_chunk, sweep_timing_events_tag_, ExePerm::EXECUTE);

        std::stringstream message_f;
        message_f << "Angleset " << angleset->GetID() << " finished on location "
                  << Chi::mpi.location_id;

        auto ev_info_f = std::make_shared<chi::ChiLog::EventInfo>(message_f.str());

        Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::SINGLE_OCCURRENCE, ev_info_f);

        scheduled_angleset++; // Schedule the next angleset
      }

      if (status != Status::FINISHED) finished = false;
    } // for each angleset rule
  }   // while not finished

  //================================================== Receive delayed data
  Chi::mpi.Barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.AngleSets())
      {
        if (angle_set->FlushSendBuffers() == Status::MESSAGES_PENDING)
          received_delayed_data = false;

        if (not angle_set->ReceiveDelayedData()) received_delayed_data = false;
      }
  }

  //================================================== Reset all
  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.AngleSets())
      angle_set->ResetSweepBuffers();

  for (auto& [bid, bndry] : angle_agg_.sim_boundaries)
  {
    if (bndry->Type() == BoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<BoundaryReflecting>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::EVENT_END);
}

void
SweepScheduler::ScheduleAlgoFIFO(SweepChunk& sweep_chunk)
{
  typedef AngleSetStatus Status;

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::EVENT_BEGIN);

  auto ev_info_i = std::make_shared<chi::ChiLog::EventInfo>(std::string("Sweep initiated"));

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::SINGLE_OCCURRENCE, ev_info_i);

  //================================================== Loop over AngleSetGroups
  AngleSetStatus completion_status = AngleSetStatus::NOT_FINISHED;
  while (completion_status == AngleSetStatus::NOT_FINISHED)
  {
    completion_status = AngleSetStatus::FINISHED;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.AngleSets())
      {
        const auto angle_set_status = angle_set->AngleSetAdvance(
          sweep_chunk, sweep_timing_events_tag_, ExecutionPermission::EXECUTE);
        if (angle_set_status == AngleSetStatus::NOT_FINISHED)
          completion_status = AngleSetStatus::NOT_FINISHED;
      } // for angleset
  }     // while not finished

  //================================================== Receive delayed data
  Chi::mpi.Barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.AngleSets())
      {
        if (angle_set->FlushSendBuffers() == Status::MESSAGES_PENDING)
          received_delayed_data = false;

        if (not angle_set->ReceiveDelayedData()) received_delayed_data = false;
      }
  }

  //================================================== Reset all
  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.AngleSets())
      angle_set->ResetSweepBuffers();

  for (auto& [bid, bndry] : angle_agg_.sim_boundaries)
  {
    if (bndry->Type() == BoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<BoundaryReflecting>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::EVENT_END);
}

void
SweepScheduler::Sweep()
{
  if (scheduler_type_ == SchedulingAlgorithm::FIRST_IN_FIRST_OUT) ScheduleAlgoFIFO(sweep_chunk_);
  else if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    ScheduleAlgoDOG(sweep_chunk_);
}

double
SweepScheduler::GetAverageSweepTime() const
{
  return Chi::log.ProcessEvent(sweep_event_tag_, chi::ChiLog::EventOperation::AVERAGE_DURATION);
}

std::vector<double>
SweepScheduler::GetAngleSetTimings()
{
  std::vector<double> info;

  double total_sweep_time =
    Chi::log.ProcessEvent(sweep_event_tag_, chi::ChiLog::EventOperation::TOTAL_DURATION);

  double total_chunk_time =
    Chi::log.ProcessEvent(sweep_timing_events_tag_[0], chi::ChiLog::EventOperation::TOTAL_DURATION);

  double ratio_sweep_to_chunk = total_chunk_time / total_sweep_time;

  info.push_back(total_sweep_time);
  info.push_back(total_chunk_time);
  info.push_back(ratio_sweep_to_chunk);

  return info;
}

void
SweepScheduler::SetDestinationPhi(std::vector<double>& in_destination_phi)
{
  sweep_chunk_.SetDestinationPhi(in_destination_phi);
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
SweepScheduler::SetDestinationPsi(std::vector<double>& in_destination_psi)
{
  sweep_chunk_.SetDestinationPsi(in_destination_psi);
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
  ZeroDestinationPsi();
  ZeroDestinationPhi();
  ZeroOutgoingDelayedPsi();
}

void
SweepScheduler::SetBoundarySourceActiveFlag(bool flag_value)
{
  sweep_chunk_.SetBoundarySourceActiveFlag(flag_value);
}

} // namespace chi_mesh::sweep_management
