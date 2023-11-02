#include "AAH_AsynComm.h"

#include "mesh/SweepUtilities/AngleSet/AngleSet.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/FLUDS/AAH_FLUDS.h"

// ###################################################################
/**Initializes delayed upstream data. This method gets called
 * when a sweep scheduler is constructed.*/
void
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::InitializeDelayedUpstreamData()
{
  const auto& spds = fluds_.GetSPDS();

  const auto num_loc_deps = spds.GetDelayedLocationDependencies().size();

  fluds_.AllocateDelayedPrelocIOutgoingPsi(num_groups_, num_angles_, num_loc_deps);

  fluds_.AllocateDelayedLocalPsi(num_groups_, num_angles_);
}