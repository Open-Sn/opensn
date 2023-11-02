#include "sweepscheduler.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**This is the entry point for sweeping.*/
void
chi_mesh::sweep_management::SweepScheduler::Sweep()
{
  if (scheduler_type_ == SchedulingAlgorithm::FIRST_IN_FIRST_OUT) ScheduleAlgoFIFO(sweep_chunk_);
  else if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    ScheduleAlgoDOG(sweep_chunk_);
}

//###################################################################
/**Get average sweep time from logging system.*/
double
chi_mesh::sweep_management::SweepScheduler::GetAverageSweepTime() const
{
  return Chi::log.ProcessEvent(sweep_event_tag_, chi::ChiLog::EventOperation::AVERAGE_DURATION);
}

//###################################################################
/**Get relevant sweep timing information.
 *
 * [0] Total sweep time
 * [1] Total chunk time
 * [2] Total chunk time / total sweep time
 * */
std::vector<double>
chi_mesh::sweep_management::SweepScheduler::GetAngleSetTimings()
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