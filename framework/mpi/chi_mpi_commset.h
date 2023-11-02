#ifndef CHITECH_CHI_MPI_COMMSET_H
#define CHITECH_CHI_MPI_COMMSET_H

#include <mpi.h>
#include "mesh/chi_mesh.h"
#include "chi_runtime.h"

namespace chi
{

//################################################################### Class def
/**Simple implementation a communicator set.
 * Definitions:
 * P = total amount of processors.
 * locI = process I in [0,P]*/
class ChiMPICommunicatorSet
{
private:
  /**A list of communicators, size P, contains a communicator for
   * only communicating with the neighbors of locI.*/
  std::vector<MPI_Comm> communicators_;
  /**A list of groupings, size P, allows mapping of the rank of locJ
   * relative to the local communicator.*/
  std::vector<MPI_Group> location_groups_;
  /**Used to translate ranks.*/
  MPI_Group world_group_;

public:
  ChiMPICommunicatorSet(std::vector<MPI_Comm>& communicators,
                        std::vector<MPI_Group>& location_groups,
                        MPI_Group& world_group)
    : communicators_(communicators), location_groups_(location_groups), world_group_(world_group)
  {
  }

  MPI_Comm LocICommunicator(int locI) const { return communicators_[locI]; }

  int MapIonJ(int locI, int locJ) const
  {
    int group_rank;
    MPI_Group_translate_ranks(world_group_, 1, &locI, location_groups_[locJ], &group_rank);

    return group_rank;
  }
};
} // namespace chi

#endif // CHITECH_CHI_MPI_COMMSET_H
