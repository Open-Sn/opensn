// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace opensn
{

AAH_ASynchronousCommunicator::AAH_ASynchronousCommunicator(
  FLUDS& fluds,
  std::size_t groupset_id,
  std::size_t angle_set_id,
  unsigned int num_groups,
  std::size_t num_angles,
  int max_mpi_message_size,
  bool csda_enabled,
  const SweepCommunicator& sweep_communicator,
  std::shared_ptr<SweepSendWindow> send_window)
  : AsynchronousCommunicator(fluds, groupset_id, angle_set_id, sweep_communicator),
    num_groups_(num_groups),
    num_angles_(num_angles),
    max_mpi_message_size_(max_mpi_message_size),
    csda_enabled_(csda_enabled),
    done_sending_(false),
    data_initialized_(false),
    upstream_data_initialized_(false),
    send_window_(std::move(send_window)),
    send_unit_size_(csda_enabled ? 2 : 1)
{
  if (not send_window_)
    throw std::invalid_argument("AAH_ASynchronousCommunicator: Missing send window.");
  // CSDA psi and psiE messages use separate tag blocks.
  num_message_tag_blocks_ = csda_enabled_ ? 2 : 1;
  BuildMessageStructure();
}

bool
AAH_ASynchronousCommunicator::DoneSending() const
{
  return done_sending_;
}

void
AAH_ASynchronousCommunicator::ClearLocalAndReceiveBuffers()
{
  fluds_.ClearLocalAndReceivePsi();
}

bool
AAH_ASynchronousCommunicator::UndelayedSendsPosted() const
{
  return sends_queued_ and next_send_unit_ == send_units_.size();
}

void
AAH_ASynchronousCommunicator::RetireCompletedSends()
{
  if (posted_requests_.empty())
    return;

  if (not mpi::test_some(posted_requests_, completed_requests_) or completed_requests_.empty())
    return;

  // MPI sets completed requests to MPI_REQUEST_NULL. Remove them so later tests skip them.
  std::erase_if(posted_requests_,
                [](const mpi::Request& request)
                { return static_cast<const MPI_Request&>(request) == MPI_REQUEST_NULL; });
  // MPI_Testsome reports the completion count as an int.
  send_window_->Release(static_cast<int>(completed_requests_.size()));
}

void
AAH_ASynchronousCommunicator::PostSendUnits(const std::vector<SendUnit>& units, std::size_t& next)
{
  const auto& comm = sweep_communicator_.GetCommunicator();
  while (next < units.size() and send_window_->TryAcquire(send_unit_size_))
  {
    const auto& [successor, m] = units[next];
    const auto& [dest, size, block_pos] = deploc_msg_data_[successor][m];
    const auto& outgoing_psi = fluds_.DeplocIOutgoingPsi()[successor];
    posted_requests_.push_back(
      comm.isend(dest, GetMessageTag(m, 0), &outgoing_psi[block_pos], size));
    // The receiver takes psi only once psiE has also arrived, so the pair is posted together.
    if (csda_enabled_)
    {
      const auto& outgoing_psiE = fluds_.DeplocIOutgoingPsiE()[successor];
      posted_requests_.push_back(
        comm.isend(dest, GetMessageTag(m, 1), &outgoing_psiE[block_pos], size));
    }
    ++next;
  }
}

void
AAH_ASynchronousCommunicator::ProgressDownstreamPsi(bool post_delayed)
{
  if (done_sending_ or not sends_queued_)
    return;

  RetireCompletedSends();
  PostSendUnits(send_units_, next_send_unit_);
  if (post_delayed)
    PostSendUnits(delayed_send_units_, next_delayed_send_unit_);

  if (next_send_unit_ == send_units_.size() and
      next_delayed_send_unit_ == delayed_send_units_.size() and posted_requests_.empty())
  {
    done_sending_ = true;
    fluds_.ClearSendPsi();
  }
}

void
AAH_ASynchronousCommunicator::Reset()
{
  done_sending_ = false;
  data_initialized_ = false;
  upstream_data_initialized_ = false;
  sends_queued_ = false;
  next_send_unit_ = 0;
  next_delayed_send_unit_ = 0;
  if (not posted_requests_.empty())
    throw std::logic_error("AAH_ASynchronousCommunicator: Reset with sweep sends outstanding.");

  for (auto& rcv_flags : preloc_msg_received_)
    rcv_flags.assign(rcv_flags.size(), false);

  for (auto& rcv_flags : delayed_preloc_msg_received_)
    rcv_flags.assign(rcv_flags.size(), false);
}

void
AAH_ASynchronousCommunicator::BuildMessageStructure()
{
  const auto& spds = fluds_.GetSPDS();
  auto* aah_fluds = dynamic_cast<AAH_FLUDS*>(&fluds_);
  if (aah_fluds == nullptr)
    throw std::runtime_error("AAH_ASynchronousCommunicator does not get AAH_FLUDS.\n");
  // Predecessor locations
  SetupMessageData(
    spds.GetLocationDependencies(),
    [this, aah_fluds](std::size_t i)
    { return aah_fluds->GetPrelocIFaceDOFCount(i) * num_groups_ * num_angles_; },
    preloc_msg_data_,
    &preloc_msg_received_,
    max_num_messages_,
    sweep_communicator_,
    max_mpi_message_size_);
  // Delayed predecessor locations
  SetupMessageData(
    spds.GetDelayedLocationDependencies(),
    [this, aah_fluds](std::size_t i)
    { return aah_fluds->GetDelayedPrelocIFaceDOFCount(i) * num_groups_ * num_angles_; },
    delayed_preloc_msg_data_,
    &delayed_preloc_msg_received_,
    max_num_messages_,
    sweep_communicator_,
    max_mpi_message_size_);
  // Successor locations
  SetupMessageData(
    spds.GetLocationSuccessors(),
    [this, aah_fluds](std::size_t i)
    { return aah_fluds->GetDeplocIFaceDOFCount(i) * num_groups_ * num_angles_; },
    deploc_msg_data_,
    nullptr,
    max_num_messages_,
    sweep_communicator_,
    max_mpi_message_size_);
  // Split the messages by whether the successor takes them during or after the sweep.
  const auto& successors = spds.GetLocationSuccessors();
  const auto& delayed_successors = spds.GetDelayedLocationSuccessors();
  send_units_.clear();
  delayed_send_units_.clear();
  for (std::size_t i = 0; i < deploc_msg_data_.size(); ++i)
  {
    const bool delayed =
      std::find(delayed_successors.begin(), delayed_successors.end(), successors[i]) !=
      delayed_successors.end();
    auto& units = delayed ? delayed_send_units_ : send_units_;
    // SetupMessageData counts the messages of a location in an int.
    const auto num_messages = static_cast<int>(deploc_msg_data_[i].size());
    for (int m = 0; m < num_messages; ++m)
      units.push_back({i, m});
  }
}

void
AAH_ASynchronousCommunicator::InitializeDelayedUpstreamData()
{
  fluds_.AllocateDelayedPrelocIOutgoingPsi();
  fluds_.AllocateDelayedLocalPsi();
}

bool
AAH_ASynchronousCommunicator::ReceiveDelayedData()
{

  const auto& spds = fluds_.GetSPDS();
  const auto& comm = sweep_communicator_.GetCommunicator();
  const std::size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();

  bool all_messages_received = true;
  for (std::size_t i = 0; i < num_delayed_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[i];

    for (int m = 0; m < delayed_preloc_msg_data_[i].size(); ++m)
    {
      const auto& [source, size, block_pos] = delayed_preloc_msg_data_[i][m];
      const int psi_tag = GetMessageTag(m, 0);
      if (not delayed_preloc_msg_received_[i][m])
      {
        const int psiE_tag = csda_enabled_ ? GetMessageTag(m, 1) : 0;
        if (not comm.iprobe(source, psi_tag) or
            (csda_enabled_ and not comm.iprobe(source, psiE_tag)))
        {
          all_messages_received = false;
          continue;
        }
        const bool psi_ok =
          not comm.recv<double>(source, psi_tag, &upstream_psi[block_pos], size).error();
        bool psiE_ok = true;
        if (csda_enabled_)
        {
          auto& upstream_psiE = fluds_.DelayedPrelocIOutgoingPsiE()[i];
          psiE_ok =
            not comm.recv<double>(source, psiE_tag, &upstream_psiE[block_pos], size).error();
        }
        if (psi_ok and psiE_ok)
          delayed_preloc_msg_received_[i][m] = true;
      }
    }
  }

  return all_messages_received;
}

AngleSetStatus
AAH_ASynchronousCommunicator::ReceiveUpstreamPsi()
{

  const auto& spds = fluds_.GetSPDS();
  const auto& comm = sweep_communicator_.GetCommunicator();
  const std::size_t num_dependencies = spds.GetLocationDependencies().size();

  // Resize FLUDS non-local incoming data
  if (not upstream_data_initialized_)
  {
    fluds_.AllocatePrelocIOutgoingPsi();
    upstream_data_initialized_ = true;
  }

  bool all_messages_received = true;
  for (std::size_t i = 0; i < num_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[i];

    for (int m = 0; m < preloc_msg_data_[i].size(); ++m)
    {
      const auto& [source, size, block_pos] = preloc_msg_data_[i][m];
      const int psi_tag = GetMessageTag(m, 0);
      if (not preloc_msg_received_[i][m])
      {
        const int psiE_tag = csda_enabled_ ? GetMessageTag(m, 1) : 0;
        if (not comm.iprobe(source, psi_tag) or
            (csda_enabled_ and not comm.iprobe(source, psiE_tag)))
        {
          all_messages_received = false;
          continue;
        }
        const bool psi_ok =
          not comm.recv<double>(source, psi_tag, &upstream_psi[block_pos], size).error();
        bool psiE_ok = true;
        if (csda_enabled_)
        {
          auto& upstream_psiE = fluds_.PrelocIOutgoingPsiE()[i];
          psiE_ok =
            not comm.recv<double>(source, psiE_tag, &upstream_psiE[block_pos], size).error();
        }
        if (psi_ok and psiE_ok)
          preloc_msg_received_[i][m] = true;
      }
    }
  }

  // Every available message is taken before returning, so that senders get their transmit
  // resources back and no matching waits on another dependency.
  return all_messages_received ? AngleSetStatus::READY_TO_EXECUTE : AngleSetStatus::RECEIVING;
}

void
AAH_ASynchronousCommunicator::SendDownstreamPsi()
{
  sends_queued_ = true;
  ProgressDownstreamPsi(false);
}

void
AAH_ASynchronousCommunicator::InitializeLocalAndDownstreamBuffers()
{
  if (not data_initialized_)
  {
    const auto& spds = fluds_.GetSPDS();

    // Resize FLUDS local outgoing data
    fluds_.AllocateInternalLocalPsi();

    // Resize FLUDS non-local outgoing data
    fluds_.AllocateOutgoingPsi();

    data_initialized_ = true;
  }
}

} // namespace opensn
