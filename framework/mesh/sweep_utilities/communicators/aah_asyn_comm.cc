#include "framework/mesh/sweep_utilities/communicators/aah_asyn_comm.h"
#include "framework/mesh/sweep_utilities/angle_set/angle_set.h"
#include "framework/mesh/sweep_utilities/spds/spds.h"
#include "framework/mesh/sweep_utilities/fluds/aah_fluds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"
#include "framework/memory_usage.h"

namespace opensn
{

AAH_ASynchronousCommunicator::AAH_ASynchronousCommunicator(FLUDS& fluds,
                                                           size_t num_groups,
                                                           size_t num_angles,
                                                           int sweep_eager_limit,
                                                           const MPICommunicatorSet& in_comm_set)
  : AsynchronousCommunicator(fluds, in_comm_set), num_groups_(num_groups), num_angles_(num_angles)
{
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;
  EAGER_LIMIT = sweep_eager_limit;

  max_num_mess = 0;

  this->BuildMessageStructure();
}

bool
AAH_ASynchronousCommunicator::DoneSending() const
{
  return done_sending;
}

void
AAH_ASynchronousCommunicator::ClearLocalAndReceiveBuffers()
{
  fluds_.ClearLocalAndReceivePsi();
}

void
AAH_ASynchronousCommunicator::ClearDownstreamBuffers()
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

void
AAH_ASynchronousCommunicator::Reset()
{
  done_sending = false;
  data_initialized = false;
  upstream_data_initialized = false;

  for (auto& message_flags : prelocI_message_received)
    message_flags.assign(message_flags.size(), false);

  for (auto& message_flags : delayed_prelocI_message_received)
    message_flags.assign(message_flags.size(), false);
}

void
AAH_ASynchronousCommunicator::BuildMessageStructure()
{
  const auto& spds = fluds_.GetSPDS();
  auto& aah_fluds = dynamic_cast<AAH_FLUDS&>(fluds_);

  // Predecessor locations
  size_t num_dependencies = spds.GetLocationDependencies().size();

  prelocI_message_count.resize(num_dependencies, 0);
  prelocI_message_size.resize(num_dependencies);
  prelocI_message_blockpos.resize(num_dependencies);
  prelocI_message_received.clear();

  for (int prelocI = 0; prelocI < num_dependencies; prelocI++)
  {
    u_ll_int num_unknowns = aah_fluds.GetPrelocIFaceDOFCount(prelocI) * num_groups_ * num_angles_;

    u_ll_int message_size;
    int message_count;
    if ((num_unknowns * 8) <= EAGER_LIMIT)
    {
      message_count = static_cast<int>(num_angles_);
      message_size = ceil((double)num_unknowns / (double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns * 8 / (double)(double)EAGER_LIMIT);
      message_size = ceil((double)num_unknowns / (double)message_count);
    }

    prelocI_message_count[prelocI] = message_count;

    u_ll_int pre_block_pos = 0;
    for (int m = 0; m < (message_count - 1); m++)
    {
      prelocI_message_size[prelocI].push_back(message_size);
      prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns > 0)
    {
      prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      prelocI_message_size[prelocI].push_back(num_unknowns);
    }

    prelocI_message_received.emplace_back(message_count, false);
  } // for prelocI

  // Delayed Predecessor locations
  size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();

  delayed_prelocI_message_count.resize(num_delayed_dependencies, 0);
  delayed_prelocI_message_size.resize(num_delayed_dependencies);
  delayed_prelocI_message_blockpos.resize(num_delayed_dependencies);
  delayed_prelocI_message_received.clear();

  for (int prelocI = 0; prelocI < num_delayed_dependencies; prelocI++)
  {
    u_ll_int num_unknowns =
      aah_fluds.GetDelayedPrelocIFaceDOFCount(prelocI) * num_groups_ * num_angles_;

    u_ll_int message_size;
    int message_count;
    if ((num_unknowns * 8) <= EAGER_LIMIT)
    {
      message_count = static_cast<int>(num_angles_);
      message_size = ceil((double)num_unknowns / (double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns * 8 / (double)(double)EAGER_LIMIT);
      message_size = ceil((double)num_unknowns / (double)message_count);
    }

    delayed_prelocI_message_count[prelocI] = message_count;

    u_ll_int pre_block_pos = 0;
    for (int m = 0; m < (message_count - 1); m++)
    {
      delayed_prelocI_message_size[prelocI].push_back(message_size);
      delayed_prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns > 0)
    {
      delayed_prelocI_message_blockpos[prelocI].push_back(pre_block_pos);
      delayed_prelocI_message_size[prelocI].push_back(num_unknowns);
    }

    delayed_prelocI_message_received.emplace_back(message_count, false);
  }

  // Successor locations
  size_t num_successors = spds.GetLocationSuccessors().size();

  deplocI_message_count.resize(num_successors, 0);
  deplocI_message_size.resize(num_successors);
  deplocI_message_blockpos.resize(num_successors);

  deplocI_message_request.clear();

  for (int deplocI = 0; deplocI < num_successors; deplocI++)
  {
    u_ll_int num_unknowns = aah_fluds.GetDeplocIFaceDOFCount(deplocI) * num_groups_ * num_angles_;

    u_ll_int message_size;
    int message_count;
    if ((num_unknowns * 8) <= EAGER_LIMIT)
    {
      message_count = static_cast<int>(num_angles_);
      message_size = ceil((double)num_unknowns / (double)message_count);
    }
    else
    {
      message_count = ceil((double)num_unknowns * 8 / (double)(double)EAGER_LIMIT);
      message_size = ceil((double)num_unknowns / (double)message_count);
    }

    deplocI_message_count[deplocI] = message_count;

    u_ll_int dep_block_pos = 0;
    for (int m = 0; m < (message_count - 1); m++)
    {
      deplocI_message_size[deplocI].push_back(message_size);
      deplocI_message_blockpos[deplocI].push_back(dep_block_pos);
      num_unknowns -= message_size;
      dep_block_pos += message_size;
    }
    if (num_unknowns > 0)
    {
      deplocI_message_blockpos[deplocI].push_back(dep_block_pos);
      deplocI_message_size[deplocI].push_back(num_unknowns);
    }

    deplocI_message_request.emplace_back(message_count, MPI_Request());
  }

  // All reduce to get maximum message count
  int angset_max_message_count = 0;
  for (size_t prelocI = 0; prelocI < num_dependencies; prelocI++)
    angset_max_message_count = std::max(prelocI_message_count[prelocI], angset_max_message_count);

  for (size_t prelocI = 0; prelocI < num_delayed_dependencies; prelocI++)
    angset_max_message_count =
      std::max(delayed_prelocI_message_count[prelocI], angset_max_message_count);

  for (size_t deplocI = 0; deplocI < num_successors; deplocI++)
    angset_max_message_count = std::max(deplocI_message_count[deplocI], angset_max_message_count);

  // Temporarily assign max_num_mess tot he local maximum
  max_num_mess = angset_max_message_count;
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

  const auto& delayed_location_dependencies = spds.GetDelayedLocationDependencies();
  const size_t num_delayed_loc_deps = delayed_location_dependencies.size();

  // Receive delayed data
  bool all_messages_received = true;
  for (size_t prelocI = 0; prelocI < num_delayed_loc_deps; prelocI++)
  {
    int locJ = delayed_location_dependencies[prelocI];

    int num_mess = delayed_prelocI_message_count[prelocI];
    for (int m = 0; m < num_mess; m++)
    {
      if (not delayed_prelocI_message_received[prelocI][m])
      {
        int message_available = 0;
        MPI_Iprobe(comm_set_.MapIonJ(locJ, opensn::mpi.location_id),
                   max_num_mess * angle_set_num + m,
                   comm_set_.LocICommunicator(opensn::mpi.location_id),
                   &message_available,
                   MPI_STATUS_IGNORE);

        if (not message_available)
        {
          all_messages_received = false;
          continue;
        }

        // Receive upstream data
        auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[prelocI];

        u_ll_int block_addr = delayed_prelocI_message_blockpos[prelocI][m];
        u_ll_int message_size = delayed_prelocI_message_size[prelocI][m];

        int error_code = MPI_Recv(&upstream_psi[block_addr],
                                  static_cast<int>(message_size),
                                  MPI_DOUBLE,
                                  comm_set_.MapIonJ(locJ, opensn::mpi.location_id),
                                  max_num_mess * angle_set_num + m,
                                  comm_set_.LocICommunicator(opensn::mpi.location_id),
                                  MPI_STATUS_IGNORE);

        delayed_prelocI_message_received[prelocI][m] = true;

        if (error_code != MPI_SUCCESS)
        {
          std::stringstream err_stream;
          err_stream << "################# Delayed receive error."
                     << " message size=" << message_size << " as_num=" << angle_set_num
                     << " num_mess=" << num_mess << " m=" << m << " error="
                     << " size="
                     << "\n";
          char error_string[BUFSIZ];
          int length_of_error_string, error_class;
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          MPI_Error_string(error_code, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          log.LogAllWarning() << err_stream.str();
        }
      } // if not message already received
    }   // for message
  }     // for delayed predecessor

  if (not all_messages_received) return false;

  return true;
}

AngleSetStatus
AAH_ASynchronousCommunicator::ReceiveUpstreamPsi(int angle_set_num)
{
  const auto& spds = fluds_.GetSPDS();

  // Resize FLUDS non-local incoming Data
  const size_t num_loc_deps = spds.GetLocationDependencies().size();
  if (!upstream_data_initialized)
  {
    fluds_.AllocatePrelocIOutgoingPsi(num_groups_, num_angles_, num_loc_deps);

    upstream_data_initialized = true;
  }

  // Assume all data is available and now try to receive all of it
  bool ready_to_execute = true;
  for (size_t prelocI = 0; prelocI < num_loc_deps; prelocI++)
  {
    int locJ = spds.GetLocationDependencies()[prelocI];

    size_t num_mess = prelocI_message_count[prelocI];
    for (int m = 0; m < num_mess; m++)
    {
      if (!prelocI_message_received[prelocI][m])
      {
        int message_available = 0;
        MPI_Iprobe(comm_set_.MapIonJ(locJ, opensn::mpi.location_id),
                   max_num_mess * angle_set_num + m,
                   comm_set_.LocICommunicator(opensn::mpi.location_id),
                   &message_available,
                   MPI_STATUS_IGNORE);

        if (not message_available)
        {
          ready_to_execute = false;
          continue;
        } // if message is not available

        // Receive upstream data
        auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[prelocI];

        u_ll_int block_addr = prelocI_message_blockpos[prelocI][m];
        u_ll_int message_size = prelocI_message_size[prelocI][m];

        int error_code = MPI_Recv(&upstream_psi[block_addr],
                                  static_cast<int>(message_size),
                                  MPI_DOUBLE,
                                  comm_set_.MapIonJ(locJ, opensn::mpi.location_id),
                                  max_num_mess * angle_set_num + m,
                                  comm_set_.LocICommunicator(opensn::mpi.location_id),
                                  MPI_STATUS_IGNORE);

        prelocI_message_received[prelocI][m] = true;

        if (error_code != MPI_SUCCESS)
        {
          std::stringstream err_stream;
          err_stream << "################# Delayed receive error."
                     << " message size=" << message_size << " as_num=" << angle_set_num
                     << " num_mess=" << num_mess << " m=" << m << " error="
                     << " size=\n";
          char error_string[BUFSIZ];
          int length_of_error_string, error_class;
          MPI_Error_class(error_code, &error_class);
          MPI_Error_string(error_class, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          MPI_Error_string(error_code, error_string, &length_of_error_string);
          err_stream << error_string << "\n";
          log.LogAllWarning() << err_stream.str();
        }
      } // if not message already received
    }   // for message

    if (!ready_to_execute) break;
  } // for predecessor

  if (!ready_to_execute) return AngleSetStatus::RECEIVING;
  else
    return AngleSetStatus::READY_TO_EXECUTE;
}

void
AAH_ASynchronousCommunicator::SendDownstreamPsi(int angle_set_num)
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
                max_num_mess * angle_set_num + m,
                comm_set_.LocICommunicator(locJ),
                &deplocI_message_request[deplocI][m]);
    } // for message
  }   // for deplocI
}

void
AAH_ASynchronousCommunicator::InitializeLocalAndDownstreamBuffers()
{
  if (!data_initialized)
  {
    const auto& spds = fluds_.GetSPDS();

    // Resize FLUDS local outgoing Data
    fluds_.AllocateInternalLocalPsi(num_groups_, num_angles_);

    // Resize FLUDS non-local outgoing Data
    const size_t num_loc_sucs = spds.GetLocationSuccessors().size();
    fluds_.AllocateOutgoingPsi(num_groups_, num_angles_, num_loc_sucs);

    // Make a memory query
    double memory_mb = GetMemoryUsageInMB();

    std::shared_ptr<Logger::EventInfo> memory_event_info =
      std::make_shared<Logger::EventInfo>(memory_mb);

    log.LogEvent(
      Logger::StdTags::MAX_MEMORY_USAGE, Logger::EventType::SINGLE_OCCURRENCE, memory_event_info);

    data_initialized = true;
  }
}

} // namespace opensn
