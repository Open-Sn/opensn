// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mpi/mpi_utils.h"
#include <cassert>

namespace opensn
{

MatchedMessage::MatchedMessage(MPI_Message message, const MPI_Status& status)
  : message_(message), status_(status)
{
}

MatchedMessage::MatchedMessage(MatchedMessage&& other) noexcept
  : message_(std::exchange(other.message_, MPI_MESSAGE_NULL)), status_(other.status_)
{
}

MatchedMessage&
MatchedMessage::operator=(MatchedMessage&& other) noexcept
{
  if (this != &other)
  {
    assert(message_ == MPI_MESSAGE_NULL);
    message_ = std::exchange(other.message_, MPI_MESSAGE_NULL);
    status_ = other.status_;
  }
  return *this;
}

int
MatchedMessage::source() const noexcept
{
  return status_.MPI_SOURCE;
}

int
MatchedMessage::tag() const noexcept
{
  return status_.MPI_TAG;
}

MatchedMessage
IProbeMatchedMessage(int source, int tag, const mpi::Communicator& comm)
{
  int message_available = 0;
  MPI_Message message = MPI_MESSAGE_NULL;
  MPI_Status status = {};
  [[maybe_unused]] const auto error_code =
    MPI_Improbe(source, tag, static_cast<MPI_Comm>(comm), &message_available, &message, &status);
  assert(error_code == MPI_SUCCESS);
  if (message_available == 0)
    return {};

  return {message, status};
}

std::vector<uint64_t>
BuildLocationExtents(uint64_t local_size, const mpi::Communicator& comm)
{
  const int process_count = comm.size();
  // Get the local vector sizes per process
  std::vector<uint64_t> local_sizes;
  comm.all_gather(local_size, local_sizes);

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  std::vector<uint64_t> extents(process_count + 1, 0);
  for (int locJ = 1; locJ < process_count; ++locJ)
    extents[locJ] = extents[locJ - 1] + local_sizes[locJ - 1];
  extents[process_count] = extents[process_count - 1] + local_sizes.back();

  return extents;
}

} // namespace opensn
