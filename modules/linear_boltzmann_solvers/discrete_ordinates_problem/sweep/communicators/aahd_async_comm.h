// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_message_struct.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include "caribou/main.hpp"
#include <vector>

namespace crb = caribou;
namespace mpi = mpicpp_lite;

namespace opensn
{

/// Interprocess communicator for AAH sweep on GPU devices.
class AAHD_ASynchronousCommunicator : public AsynchronousCommunicator
{
public:
  AAHD_ASynchronousCommunicator(FLUDS& fluds,
                                std::size_t num_groups,
                                std::size_t num_angles,
                                int max_mpi_message_size,
                                const MPICommunicatorSet& comm_set);

  int GetMaxNumMessages() const { return max_num_messages_; }

  void SetMaxNumMessages(int count) { max_num_messages_ = count; }

  /**
   * Initialize delayed upstream data.
   * This method gets called when a sweep scheduler is constructed.
   */
  void InitializeDelayedUpstreamData();

  /// Allocate local and non-local outgoing memory.
  void InitializeLocalAndDownstreamBuffers(crb::Stream& stream);

  /// Block until all upstream dependencies have been met.
  void ReceiveUpstreamPsi(int angle_set_num, crb::Stream& stream);

  /// Send non-local outgoing psi.
  void SendDownstreamPsi(int angle_set_num, crb::Stream& stream);

  /// Receive delayed data from successor locations.
  void ReceiveDelayedData(int angle_set_num);

  /// Wait until all outgoing messages have been sent and all delayed incoming messages have been received.
  void Wait();

  /// Deallocate memory for local and non-local incoming psi after the sweep.
  void ClearLocalAndReceiveBuffers(crb::Stream& stream);

  /// Clear flags in preperation for another sweep.
  void Reset();

protected:
  /**
   * Build message structure.
   * Message structure for tracking incoming, delayed incoming and outgoing non-local fluxes.
   */
  void BuildMessageStructure();

private:
  int max_num_messages_;
  int max_mpi_message_size_;

  std::vector<mpi::Request> preloc_msg_request_;
  std::vector<std::vector<AAH_MessageDetails>> preloc_msg_data_;

  std::vector<mpi::Request> delayed_preloc_msg_request_;
  std::vector<std::vector<AAH_MessageDetails>> delayed_preloc_msg_data_;

  std::vector<mpi::Request> deploc_msg_request_;
  std::vector<std::vector<AAH_MessageDetails>> deploc_msg_data_;
};

} // namespace opensn
