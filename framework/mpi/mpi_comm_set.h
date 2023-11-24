#pragma once

#include "framework/mesh/mesh.h"
#include "framework/runtime.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

namespace opensn
{

/**Simple implementation a communicator set.
 * Definitions:
 * P = total amount of processors.
 * locI = process I in [0,P]*/
class MPICommunicatorSet
{
private:
  /**A list of communicators, size P, contains a communicator for
   * only communicating with the neighbors of locI.*/
  std::vector<mpi::Communicator> communicators_;
  /**A list of groupings, size P, allows mapping of the rank of locJ
   * relative to the local communicator.*/
  std::vector<mpi::Group> location_groups_;
  /**Used to translate ranks.*/
  mpi::Group world_group_;

public:
  MPICommunicatorSet(std::vector<mpi::Communicator>& communicators,
                     std::vector<mpi::Group>& location_groups,
                     mpi::Group& world_group)
    : communicators_(communicators), location_groups_(location_groups), world_group_(world_group)
  {
  }

  const mpi::Communicator& LocICommunicator(int locI) const { return communicators_[locI]; }

  int MapIonJ(int locI, int locJ) const
  {
    return world_group_.translate_rank(locI, location_groups_[locJ]);
  }
};

} // namespace opensn
