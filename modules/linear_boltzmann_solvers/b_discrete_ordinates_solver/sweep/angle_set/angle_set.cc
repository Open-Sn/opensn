#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/angle_set/angle_set.h"
#include "framework/logging/log_exceptions.h"

namespace opensn
{
namespace lbs
{

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

size_t
AngleSet::GetID() const
{
  return id_;
}

const SPDS&
AngleSet::GetSPDS() const
{
  return spds_;
}

FLUDS&
AngleSet::GetFLUDS()
{
  return *fluds_;
}

size_t
AngleSet::GetRefGroupSubset() const
{
  return ref_group_subset_;
}

const std::vector<size_t>&
AngleSet::GetAngleIndices() const
{
  return angles_;
}

std::map<uint64_t, std::shared_ptr<SweepBndry>>&
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

} // namespace lbs
} // namespace opensn
