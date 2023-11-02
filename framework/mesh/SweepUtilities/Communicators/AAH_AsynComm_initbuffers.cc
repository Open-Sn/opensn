#include "AAH_AsynComm.h"

#include "mesh/SweepUtilities/AngleSet/AngleSet.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/AAH_FLUDS.h"

#include "chi_runtime.h"
#include "console/chi_console.h"
#include "chi_log.h"

//###################################################################
/** This is the final level of initialization before a sweep-chunk executes.
 * Once all upstream dependencies are met and if the sweep scheduler places
 * this angleset as "ready-to-execute", then the angle-set will call this
 * method. It is also fairly important in terms of memory to only allocate
 * these chunks of memory when actually ready to use them since they form the
 * majority of memory usage.*/
void
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::InitializeLocalAndDownstreamBuffers()
{
  if (!data_initialized)
  {
    const auto& spds = fluds_.GetSPDS();

    //============================ Resize FLUDS local outgoing Data
    fluds_.AllocateInternalLocalPsi(num_groups_, num_angles_);

    //============================ Resize FLUDS non-local outgoing Data
    const size_t num_loc_sucs = spds.GetLocationSuccessors().size();
    fluds_.AllocateOutgoingPsi(num_groups_, num_angles_, num_loc_sucs);

    //================================================ Make a memory query
    double memory_mb = chi::Console::GetMemoryUsageInMB();

    std::shared_ptr<chi::ChiLog::EventInfo> memory_event_info =
      std::make_shared<chi::ChiLog::EventInfo>(memory_mb);

    Chi::log.LogEvent(chi::ChiLog::StdTags::MAX_MEMORY_USAGE,
                      chi::ChiLog::EventType::SINGLE_OCCURRENCE,
                      memory_event_info);

    data_initialized = true;
  }
}