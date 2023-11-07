#include "framework/mesh/SweepUtilities/AngleSet/AngleSet.h"

#include "framework/logging/chi_log_exceptions.h"

namespace chi_mesh::sweep_management
{

/**AngleSet constructor.*/
AngleSet::AngleSet(size_t id,
                   size_t num_groups,
                   const SPDS& spds,
                   std::shared_ptr<FLUDS>& fluds,
                   const std::vector<size_t>& angle_indices,
                   std::map<uint64_t, std::shared_ptr<SweepBndry>>& sim_boundaries,
                   const size_t in_ref_subset)
  : id_(id),
    num_grps(num_groups),
    spds_(spds),
    fluds_(fluds),
    angles_(angle_indices),
    ref_boundaries_(sim_boundaries),
    ref_group_subset_(in_ref_subset)
{
}

/**Returns the angleset's unique id.*/
size_t
AngleSet::GetID() const
{
  return id_;
}

/**Returns a reference to the associated spds.*/
const SPDS&
AngleSet::GetSPDS() const
{
  return spds_;
}

/**Returns a reference to the associated fluds_.*/
FLUDS&
AngleSet::GetFLUDS()
{
  return *fluds_;
}

/**Return the reference group subset number.*/
size_t
AngleSet::GetRefGroupSubset() const
{
  return ref_group_subset_;
}

/**Returns the angle indices associated with this angleset.*/
const std::vector<size_t>&
AngleSet::GetAngleIndices() const
{
  return angles_;
}

/**Returns the angle indices associated with this angleset.*/
std::map<uint64_t, AngleSet::SweepBndryPtr>&
AngleSet::GetBoundaries()
{
  return ref_boundaries_;
}

size_t
AngleSet::GetNumGroups() const
{
  return num_grps;
}
size_t
AngleSet::GetNumAngles() const
{
  return angles_.size();
}

AsynchronousCommunicator*
AngleSet::GetCommunicator()
{
  ChiLogicalError("Method not implemented");
}

} // namespace chi_mesh::sweep_management
