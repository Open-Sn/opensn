#include "sweepscheduler.h"

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Applies a First-In-First-Out sweep scheduling.*/
void chi_mesh::sweep_management::SweepScheduler::ScheduleAlgoFIFO(
  SweepChunk& sweep_chunk)
{
  typedef AngleSetStatus Status;

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::EVENT_BEGIN);

  auto ev_info_i =
    std::make_shared<chi::ChiLog::EventInfo>(std::string("Sweep initiated"));

  Chi::log.LogEvent(
    sweep_event_tag_, chi::ChiLog::EventType::SINGLE_OCCURRENCE, ev_info_i);

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
      }// for angleset
  }// while not finished

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

        if (not angle_set->ReceiveDelayedData())
          received_delayed_data = false;
      }
  }

  //================================================== Reset all
  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.AngleSets())
      angle_set->ResetSweepBuffers();

  for (auto& [bid, bndry] : angle_agg_.sim_boundaries)
  {
    if (bndry->Type() == chi_mesh::sweep_management::BoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<
        chi_mesh::sweep_management::BoundaryReflecting>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }

  Chi::log.LogEvent(sweep_event_tag_, chi::ChiLog::EventType::EVENT_END);
}