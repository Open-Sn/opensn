#pragma once

#include "opensn/framework/mesh/SweepUtilities/sweep_namespace.h"
#include "opensn/framework/mpi/chi_mpi.h"

#include "opensn/framework/mesh/SweepUtilities/Communicators/AsyncComm.h"

typedef unsigned long long int u_ll_int;

namespace chi
{
class ChiMPICommunicatorSet;
}

namespace chi_mesh::sweep_management
{

class FLUDS;

/**
 * Handles the swift communication of interprocess communication related to sweeping.
 */
class AAH_ASynchronousCommunicator : public AsynchronousCommunicator
{
private:
  const size_t num_groups_;
  const size_t num_angles_;

  bool done_sending;
  bool data_initialized;
  bool upstream_data_initialized;

  u_ll_int EAGER_LIMIT = 32000;

  std::vector<int> prelocI_message_count;
  std::vector<int> deplocI_message_count;
  std::vector<int> delayed_prelocI_message_count;

  std::vector<std::vector<u_ll_int>> prelocI_message_size;
  std::vector<std::vector<u_ll_int>> deplocI_message_size;
  std::vector<std::vector<u_ll_int>> delayed_prelocI_message_size;

  std::vector<std::vector<u_ll_int>> prelocI_message_blockpos;
  std::vector<std::vector<u_ll_int>> deplocI_message_blockpos;
  std::vector<std::vector<u_ll_int>> delayed_prelocI_message_blockpos;

  std::vector<std::vector<bool>> prelocI_message_received;
  std::vector<std::vector<bool>> delayed_prelocI_message_received;

  std::vector<std::vector<MPI_Request>> deplocI_message_request;

public:
  int max_num_mess;

  AAH_ASynchronousCommunicator(FLUDS& fluds,
                               size_t num_groups,
                               size_t num_angles,
                               int sweep_eager_limit,
                               const chi::ChiMPICommunicatorSet& in_comm_set);
  /**
   * Returns the private flag done_sending.
   */
  bool DoneSending() const;

  /**
   * Initializes delayed upstream data. This method gets called
   * when a sweep scheduler is constructed.
   */
  void InitializeDelayedUpstreamData();

  /**
   * This is the final level of initialization before a sweep-chunk executes.
   * Once all upstream dependencies are met and if the sweep scheduler places
   * this angleset as "ready-to-execute", then the angle-set will call this
   * method. It is also fairly important in terms of memory to only allocate
   * these chunks of memory when actually ready to use them since they form the
   * majority of memory usage.
   */
  void InitializeLocalAndDownstreamBuffers();

  /**
   * Sends downstream psi. This method gets called after a sweep chunk has executed
   */
  void SendDownstreamPsi(int angle_set_num);

  /**
   * Receives delayed data from successor locations.
   */
  bool ReceiveDelayedData(int angle_set_num);

  /**
   * Sends downstream psi.
   */
  void ClearDownstreamBuffers();

  /**
   * Check if all upstream dependencies have been met and receives
   * it as it becomes available.
   */
  AngleSetStatus ReceiveUpstreamPsi(int angle_set_num);

  /**
   * Receive all upstream Psi. This method is called from within
   * an advancement of an angleset, right after execution.
   */
  void ClearLocalAndReceiveBuffers();

  /**
   * Clear flags in preperation for another sweep.
   */
  void Reset();

protected:
  /**
   * Builds message structure.
   *
   * Outgoing and incoming data needs to be sub-divided into messages
   * each of which is smaller than the MPI eager-limit. There are
   * three parts to this: predecessors, delayed-predecessors and successors.
   *
   * This method gets called by an angleset that subscribes to this
   * sweepbuffer.
   */
  void BuildMessageStructure();
};

} // namespace chi_mesh::sweep_management
