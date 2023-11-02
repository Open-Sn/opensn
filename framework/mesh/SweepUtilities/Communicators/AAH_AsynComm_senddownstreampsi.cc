#include "AAH_AsynComm.h"

#include "mesh/SweepUtilities/AngleSet/AngleSet.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"

#include "mpi/chi_mpi_commset.h"

//###################################################################
/**Sends downstream psi. This method gets called after a sweep chunk has
 * executed */
void
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::SendDownstreamPsi(int angle_set_num)
{
  const auto& spds = fluds_.GetSPDS();

  const auto& location_successors = spds.GetLocationSuccessors();

  const size_t num_successors = location_successors.size();
  for (size_t deplocI = 0; deplocI < num_successors; deplocI++)
  {
    int locJ = location_successors[deplocI];

    int num_mess = deplocI_message_count[deplocI];
    for (int m = 0; m < num_mess; m++)
    {
      u_ll_int block_addr = deplocI_message_blockpos[deplocI][m];
      u_ll_int message_size = deplocI_message_size[deplocI][m];

      const auto& outgoing_psi = fluds_.DeplocIOutgoingPsi()[deplocI];

      MPI_Isend(&outgoing_psi[block_addr],
                static_cast<int>(message_size),
                MPI_DOUBLE,
                comm_set_.MapIonJ(locJ, locJ),
                max_num_mess * angle_set_num + m, // tag
                comm_set_.LocICommunicator(locJ),
                &deplocI_message_request[deplocI][m]);
    } // for message
  }   // for deplocI
}