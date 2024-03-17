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
 * Handles the swift communication of interprocess communication related to sweeping.
 */
class AAH_ASynchronousCommunicator : public AsynchronousCommunicator
{
private:
  int num_groups_;
  int num_angles_;
  int max_num_messages_;
  int max_mpi_message_size_;
  bool done_sending_;
  bool data_initialized_;
  bool upstream_data_initialized_;

  std::vector<int> prelocI_message_count_;
  std::vector<int> deplocI_message_count_;
  std::vector<int> delayed_prelocI_message_count_;

  std::vector<std::vector<size_t>> prelocI_message_size_;
  std::vector<std::vector<size_t>> deplocI_message_size_;
  std::vector<std::vector<size_t>> delayed_prelocI_message_size_;

  std::vector<std::vector<size_t>> prelocI_message_blockpos_;
  std::vector<std::vector<size_t>> deplocI_message_blockpos_;
  std::vector<std::vector<size_t>> delayed_prelocI_message_blockpos_;

  std::vector<std::vector<bool>> prelocI_message_received_;
  std::vector<std::vector<bool>> delayed_prelocI_message_received_;

  std::vector<std::vector<mpi::Request>> deplocI_message_request_;

public:
  AAH_ASynchronousCommunicator(FLUDS& fluds,
                               size_t num_groups,
                               size_t num_angles,
                               int max_mpi_message_size,
                               const MPICommunicatorSet& in_comm_set);

  int GetMaxNumMessages() const { return max_num_messages_; }

  void SetMaxNumMessages(int count) { max_num_messages_ = count; }

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
   * each of which is smaller than the maximum message size. There are
   * three parts to this: predecessors, delayed-predecessors and successors.
   *
   * This method gets called by an angleset that subscribes to this
   * sweepbuffer.
   */
  void BuildMessageStructure();
};

} // namespace lbs
} // namespace opensn
