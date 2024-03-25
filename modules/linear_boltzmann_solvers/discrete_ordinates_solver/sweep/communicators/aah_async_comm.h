#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/sweep.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

namespace opensn
{

class MPICommunicatorSet;

namespace lbs
{

class FLUDS;

/**
 * Handles interprocess communication related to sweeping.
 */
class AAH_ASynchronousCommunicator : public AsynchronousCommunicator
{
private:
  size_t num_groups_;
  size_t num_angles_;
  size_t max_num_messages_;
  size_t max_mpi_message_size_;
  bool done_sending_;
  bool data_initialized_;
  bool upstream_data_initialized_;

  std::vector<std::vector<bool>> preloc_msg_received_;
  std::vector<std::vector<std::tuple<int, size_t, size_t>>> preloc_msg_data_;

  std::vector<std::vector<bool>> delayed_preloc_msg_received_;
  std::vector<std::vector<std::tuple<int, size_t, size_t>>> delayed_preloc_msg_data_;

  std::vector<mpi::Request> deploc_msg_request_;
  std::vector<std::vector<std::tuple<int, size_t, size_t>>> deploc_msg_data_;

protected:
  /**
   * Builds message structure.
   *
   * Outgoing and incoming data needs to be sub-divided into messages
   * each of which is smaller than the maximum message size. There are
   * three parts to this: predecessors, delayed-predecessors and successors.
   *
   * This method gets called by an angleset that subscribes to this
   * sweepbuffer.
   */
  void BuildMessageStructure();

public:
  AAH_ASynchronousCommunicator(FLUDS& fluds,
                               size_t num_groups,
                               size_t num_angles,
                               size_t max_mpi_message_size,
                               const MPICommunicatorSet& comm_set);

  size_t GetMaxNumMessages() const { return max_num_messages_; }

  void SetMaxNumMessages(size_t count) { max_num_messages_ = count; }

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
};

} // namespace lbs
} // namespace opensn
