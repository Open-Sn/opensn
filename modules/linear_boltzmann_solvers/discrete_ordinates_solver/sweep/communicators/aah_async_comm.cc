#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/communicators/aah_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/aah_fluds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/logging/log.h"
#include "framework/memory_usage.h"
#include "framework/runtime.h"

namespace opensn
{
namespace lbs
{

AAH_ASynchronousCommunicator::AAH_ASynchronousCommunicator(FLUDS& fluds,
                                                           size_t num_groups,
                                                           size_t num_angles,
                                                           int max_mpi_message_size,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    num_groups_(num_groups),
    num_angles_(num_angles),
    max_num_messages_(0),
    max_mpi_message_size_(max_mpi_message_size),
    done_sending_(false),
    data_initialized_(false),
    upstream_data_initialized_(false)
{
  this->BuildMessageStructure();
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

void
AAH_ASynchronousCommunicator::ClearDownstreamBuffers()
{
  if (done_sending_)
    return;

  done_sending_ = true;
  for (auto& locI_requests : deplocI_message_request_)
  {
    for (auto& request : locI_requests)
    {
      if (not mpi::test(request))
      {
        done_sending_ = false;
        return;
      }
    }
  }

  fluds_.ClearSendPsi();
}

void
AAH_ASynchronousCommunicator::Reset()
{
  done_sending_ = false;
  data_initialized_ = false;
  upstream_data_initialized_ = false;

  for (auto& message_flags : prelocI_message_received_)
    message_flags.assign(message_flags.size(), false);

  for (auto& message_flags : delayed_prelocI_message_received_)
    message_flags.assign(message_flags.size(), false);
}

void
AAH_ASynchronousCommunicator::BuildMessageStructure()
{
  const auto& spds = fluds_.GetSPDS();
  auto& fluds = dynamic_cast<AAH_FLUDS&>(fluds_);

  // Predecessor locations
  size_t num_dependencies = spds.GetLocationDependencies().size();

  prelocI_message_count_.resize(num_dependencies, 0);
  prelocI_message_size_.resize(num_dependencies);
  prelocI_message_blockpos_.resize(num_dependencies);
  prelocI_message_received_.clear();

  for (int prelocI = 0; prelocI < num_dependencies; ++prelocI)
  {
    size_t num_unknowns = fluds.GetPrelocIFaceDOFCount(prelocI) * num_groups_ * num_angles_;
    size_t message_count = num_angles_;
    if (num_unknowns * 8 > max_mpi_message_size_)
      message_count = ((num_unknowns * 8) + (max_mpi_message_size_ - 1)) / max_mpi_message_size_;
    size_t message_size = (num_unknowns + (message_count - 1)) / message_count;

    prelocI_message_count_[prelocI] = message_count;
    prelocI_message_size_[prelocI].reserve(message_count);
    prelocI_message_blockpos_[prelocI].reserve(message_count);

    size_t pre_block_pos = 0;
    for (int m = 0; m < message_count - 1; ++m)
    {
      prelocI_message_size_[prelocI].push_back(message_size);
      prelocI_message_blockpos_[prelocI].push_back(pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns > 0)
    {
      prelocI_message_blockpos_[prelocI].push_back(pre_block_pos);
      prelocI_message_size_[prelocI].push_back(num_unknowns);
    }

    prelocI_message_received_.emplace_back(message_count, false);
  }

  // Delayed Predecessor locations
  size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();

  delayed_prelocI_message_count_.resize(num_delayed_dependencies, 0);
  delayed_prelocI_message_size_.resize(num_delayed_dependencies);
  delayed_prelocI_message_blockpos_.resize(num_delayed_dependencies);
  delayed_prelocI_message_received_.clear();

  for (int prelocI = 0; prelocI < num_delayed_dependencies; ++prelocI)
  {
    size_t num_unknowns = fluds.GetDelayedPrelocIFaceDOFCount(prelocI) * num_groups_ * num_angles_;
    size_t message_count = num_angles_;
    if (num_unknowns * 8 > max_mpi_message_size_)
      message_count = ((num_unknowns * 8) + (max_mpi_message_size_ - 1)) / max_mpi_message_size_;
    size_t message_size = (num_unknowns + (message_count - 1)) / message_count;

    delayed_prelocI_message_count_[prelocI] = message_count;
    delayed_prelocI_message_size_[prelocI].reserve(message_count);
    delayed_prelocI_message_blockpos_[prelocI].reserve(message_count);

    size_t pre_block_pos = 0;
    for (int m = 0; m < message_count - 1; ++m)
    {
      delayed_prelocI_message_size_[prelocI].push_back(message_size);
      delayed_prelocI_message_blockpos_[prelocI].push_back(pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns > 0)
    {
      delayed_prelocI_message_blockpos_[prelocI].push_back(pre_block_pos);
      delayed_prelocI_message_size_[prelocI].push_back(num_unknowns);
    }

    delayed_prelocI_message_received_.emplace_back(message_count, false);
  }

  // Successor locations
  size_t num_successors = spds.GetLocationSuccessors().size();

  deplocI_message_count_.resize(num_successors, 0);
  deplocI_message_size_.resize(num_successors);
  deplocI_message_blockpos_.resize(num_successors);
  deplocI_message_request_.clear();

  for (int deplocI = 0; deplocI < num_successors; ++deplocI)
  {
    size_t num_unknowns = fluds.GetDeplocIFaceDOFCount(deplocI) * num_groups_ * num_angles_;
    size_t message_count = num_angles_;
    if (num_unknowns * 8 > max_mpi_message_size_)
      message_count = ((num_unknowns * 8) + (max_mpi_message_size_ - 1)) / max_mpi_message_size_;
    size_t message_size = (num_unknowns + (message_count - 1)) / message_count;

    deplocI_message_count_[deplocI] = message_count;
    deplocI_message_size_[deplocI].reserve(message_count);
    deplocI_message_blockpos_[deplocI].reserve(message_count);

    size_t dep_block_pos = 0;
    for (int m = 0; m < message_count - 1; ++m)
    {
      deplocI_message_size_[deplocI].push_back(message_size);
      deplocI_message_blockpos_[deplocI].push_back(dep_block_pos);
      num_unknowns -= message_size;
      dep_block_pos += message_size;
    }
    if (num_unknowns > 0)
    {
      deplocI_message_blockpos_[deplocI].push_back(dep_block_pos);
      deplocI_message_size_[deplocI].push_back(num_unknowns);
    }

    deplocI_message_request_.emplace_back(message_count, mpi::Request());
  }

  // All reduce to get maximum message count
  int angset_max_message_count = 0;
  for (size_t prelocI = 0; prelocI < num_dependencies; ++prelocI)
    angset_max_message_count = std::max(prelocI_message_count_[prelocI], angset_max_message_count);

  for (size_t prelocI = 0; prelocI < num_delayed_dependencies; ++prelocI)
    angset_max_message_count =
      std::max(delayed_prelocI_message_count_[prelocI], angset_max_message_count);

  for (size_t deplocI = 0; deplocI < num_successors; ++deplocI)
    angset_max_message_count = std::max(deplocI_message_count_[deplocI], angset_max_message_count);

  // Temporarily assign max_num_messages_ to the local maximum
  max_num_messages_ = angset_max_message_count;
}

void
AAH_ASynchronousCommunicator::InitializeDelayedUpstreamData()
{
  const auto& spds = fluds_.GetSPDS();
  const auto num_loc_deps = spds.GetDelayedLocationDependencies().size();
  fluds_.AllocateDelayedPrelocIOutgoingPsi(num_groups_, num_angles_, num_loc_deps);
  fluds_.AllocateDelayedLocalPsi(num_groups_, num_angles_);
}

bool
AAH_ASynchronousCommunicator::ReceiveDelayedData(int angle_set_num)
{
  const auto& spds = fluds_.GetSPDS();
  auto& delayed_location_dependencies = spds.GetDelayedLocationDependencies();
  const size_t num_delayed_loc_deps = delayed_location_dependencies.size();

  // Receive delayed data
  bool all_messages_received = true;
  for (size_t prelocI = 0; prelocI < num_delayed_loc_deps; ++prelocI)
  {
    int locJ = delayed_location_dependencies[prelocI];

    for (int m = 0; m < delayed_prelocI_message_count_[prelocI]; ++m)
    {
      if (not delayed_prelocI_message_received_[prelocI][m])
      {
        auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
        auto source_rank = comm_set_.MapIonJ(locJ, opensn::mpi_comm.rank());
        auto tag = max_num_messages_ * angle_set_num + m;
        auto message_available = comm.iprobe(source_rank, tag);
        if (not message_available)
        {
          all_messages_received = false;
          continue;
        }

        // Receive upstream data
        auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[prelocI];
        size_t block_addr = delayed_prelocI_message_blockpos_[prelocI][m];
        size_t message_size = delayed_prelocI_message_size_[prelocI][m];
        comm.recv(source_rank, tag, &upstream_psi[block_addr], message_size);

        delayed_prelocI_message_received_[prelocI][m] = true;
      } // if not message already received
    }   // for message
  }     // for delayed predecessor

  if (not all_messages_received)
    return false;

  return true;
}

AngleSetStatus
AAH_ASynchronousCommunicator::ReceiveUpstreamPsi(int angle_set_num)
{
  const auto& spds = fluds_.GetSPDS();

  // Resize FLUDS non-local incoming Data
  const size_t num_loc_deps = spds.GetLocationDependencies().size();
  if (not upstream_data_initialized_)
  {
    fluds_.AllocatePrelocIOutgoingPsi(num_groups_, num_angles_, num_loc_deps);
    upstream_data_initialized_ = true;
  }

  // Assume all data is available and now try to receive all of it
  bool ready_to_execute = true;
  for (size_t prelocI = 0; prelocI < num_loc_deps; ++prelocI)
  {
    int locJ = spds.GetLocationDependencies()[prelocI];

    for (int m = 0; m < prelocI_message_count_[prelocI]; ++m)
    {
      if (not prelocI_message_received_[prelocI][m])
      {
        auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
        auto source = comm_set_.MapIonJ(locJ, opensn::mpi_comm.rank());
        auto tag = max_num_messages_ * angle_set_num + m;
        if (not comm.iprobe(source, tag))
        {
          ready_to_execute = false;
          continue;
        } // if message is not available

        // Receive upstream data
        auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[prelocI];
        size_t block_addr = prelocI_message_blockpos_[prelocI][m];
        size_t message_size = prelocI_message_size_[prelocI][m];
        comm.recv(source, tag, &upstream_psi[block_addr], message_size);

        prelocI_message_received_[prelocI][m] = true;
      } // if not message already received
    }   // for message

    if (not ready_to_execute)
      return AngleSetStatus::RECEIVING;
  } // for predecessor

  return AngleSetStatus::READY_TO_EXECUTE;
}

void
AAH_ASynchronousCommunicator::SendDownstreamPsi(int angle_set_num)
{
  const auto& spds = fluds_.GetSPDS();
  const auto& location_successors = spds.GetLocationSuccessors();
  const size_t num_successors = location_successors.size();

  for (size_t deplocI = 0; deplocI < num_successors; ++deplocI)
  {
    int locJ = location_successors[deplocI];

    for (int m = 0; m < deplocI_message_count_[deplocI]; ++m)
    {
      size_t block_addr = deplocI_message_blockpos_[deplocI][m];
      size_t message_size = deplocI_message_size_[deplocI][m];
      const auto& outgoing_psi = fluds_.DeplocIOutgoingPsi()[deplocI];
      auto& comm = comm_set_.LocICommunicator(locJ);
      auto dest = comm_set_.MapIonJ(locJ, locJ);
      auto tag = max_num_messages_ * angle_set_num + m;
      deplocI_message_request_[deplocI][m] =
        comm.isend(dest, tag, &outgoing_psi[block_addr], message_size);
    } // for message
  }   // for deplocI
}

void
AAH_ASynchronousCommunicator::InitializeLocalAndDownstreamBuffers()
{
  if (not data_initialized_)
  {
    const auto& spds = fluds_.GetSPDS();

    // Resize FLUDS local outgoing Data
    fluds_.AllocateInternalLocalPsi(num_groups_, num_angles_);

    // Resize FLUDS non-local outgoing Data
    fluds_.AllocateOutgoingPsi(num_groups_, num_angles_, spds.GetLocationSuccessors().size());

    // Make a memory query
    double memory_mb = GetMemoryUsageInMB();

    std::shared_ptr<Logger::EventInfo> memory_event_info =
      std::make_shared<Logger::EventInfo>(memory_mb);

    log.LogEvent(
      Logger::StdTags::MAX_MEMORY_USAGE, Logger::EventType::SINGLE_OCCURRENCE, memory_event_info);

    data_initialized_ = true;
  }
}

} // namespace lbs
} // namespace opensn
