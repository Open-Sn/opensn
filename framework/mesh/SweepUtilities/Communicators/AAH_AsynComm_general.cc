#include "AAH_AsynComm.h"

#include "mesh/SweepUtilities/AngleSet/AngleSet.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"

// ###################################################################
/**Constructor.*/
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::AAH_ASynchronousCommunicator(
  FLUDS& fluds,
  size_t num_groups,
  size_t num_angles,
  int sweep_eager_limit,
  const chi::ChiMPICommunicatorSet& in_comm_set)
  : AsynchronousCommunicator(fluds, in_comm_set), num_groups_(num_groups), num_angles_(num_angles)
{
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;
  EAGER_LIMIT = sweep_eager_limit;

  max_num_mess = 0;

  this->BuildMessageStructure();
}

// ###################################################################
/**Returns the private flag done_sending.*/
bool
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::DoneSending() const
{
  return done_sending;
}

// ###################################################################
/**Receive all upstream Psi. This method is called from within
 * an advancement of an angleset, right after execution.*/
void
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::ClearLocalAndReceiveBuffers()
{
  fluds_.ClearLocalAndReceivePsi();
}

// ###################################################################
/**Sends downstream psi.*/
void
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::ClearDownstreamBuffers()
{
  if (done_sending) return;

  done_sending = true;
  for (auto& locI_requests : deplocI_message_request)
    for (auto& request : locI_requests)
    {
      int message_sent;
      MPI_Test(&request, &message_sent, MPI_STATUS_IGNORE);
      if (not message_sent)
      {
        done_sending = false;
        return;
      }
    }

  if (done_sending) fluds_.ClearSendPsi();
}

// ###################################################################
/**Clear flags in preperation for another sweep.*/
void
chi_mesh::sweep_management::AAH_ASynchronousCommunicator::Reset()
{
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;

  for (auto& message_flags : prelocI_message_received)
    message_flags.assign(message_flags.size(), false);

  for (auto& message_flags : delayed_prelocI_message_received)
    message_flags.assign(message_flags.size(), false);
}