// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/runtime.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace opensn
{

/// The error message for a thrown object that is not a std::exception.
inline constexpr const char* non_standard_exception_message =
  "an exception that is not a std::exception";

/// The message of `error`, or a placeholder if it has none, so that a failure is never empty.
inline std::string
DescribeError(const std::exception& error)
{
  std::string message = error.what();
  return message.empty() ? std::string("an error without a message") : message;
}

/**
 * Run `function`, which makes no MPI calls, and return the message of anything it throws: the
 * message of a std::exception (never empty), a fixed message for any other thrown object, or an
 * empty string if it returns normally. Pass the result to ThrowIfAnyRankFailed so that a failure on
 * one rank, including std::bad_alloc, is raised on every rank.
 */
template <typename Function>
std::string
CaptureError(const Function& function)
{
  try
  {
    function();
  }
  catch (const std::exception& exception)
  {
    return DescribeError(exception);
  }
  catch (...)
  {
    return non_standard_exception_message;
  }
  return {};
}

/// The longest error message, in bytes, that FindFirstRankError passes between ranks.
inline constexpr std::size_t max_collective_error_length = 4096;

/**
 * Find the lowest rank of `comm` whose `local_error` is not empty. Returns that rank and its error
 * on every rank, or nullopt on every rank if no rank failed. Collective. The message is passed in
 * a fixed-size buffer and cut to max_collective_error_length bytes, so no rank allocates memory
 * until every rank has finished communicating.
 */
inline std::optional<std::pair<int, std::string>>
FindFirstRankError(const mpi::Communicator& comm, const std::string& local_error)
{
  const int comm_size = comm.size();
  int failure_rank = local_error.empty() ? comm_size : comm.rank();
  comm.all_reduce(failure_rank, mpi::op::min<int>());
  if (failure_rank == comm_size)
    return std::nullopt;

  std::array<char, max_collective_error_length> buffer{};
  int length = 0;
  if (comm.rank() == failure_rank)
  {
    length = static_cast<int>(std::min(local_error.size(), buffer.size()));
    std::memcpy(buffer.data(), local_error.data(), static_cast<std::size_t>(length));
  }
  comm.broadcast(length, failure_rank);
  if (length > 0)
    comm.broadcast(buffer.data(), length, failure_rank);
  return std::make_pair(failure_rank, std::string(buffer.data(), static_cast<std::size_t>(length)));
}

/**
 * Throw `Exception` on every rank of `comm` if any rank's `local_error` is not empty. The message,
 * "<operation> failed on rank <r>: <error>", names the lowest failing rank and its error on every
 * rank. Collective.
 */
template <typename Exception = std::runtime_error>
void
ThrowIfAnyRankFailed(const mpi::Communicator& comm,
                     const std::string& local_error,
                     const std::string_view operation)
{
  const auto failure = FindFirstRankError(comm, local_error);
  if (failure.has_value())
    throw Exception(std::string(operation) + " failed on rank " + std::to_string(failure->first) +
                    ": " + failure->second);
}

/**
 * Runs rank-local work that makes no MPI calls and propagates a failure to every rank. If the work
 * throws on any rank, every rank throws std::runtime_error with the message from the lowest failing
 * rank, so no rank is left waiting in a later collective. Collective.
 */
template <typename Function>
void
RunWithCollectiveErrors(const Function& function)
{
  const auto failure = FindFirstRankError(mpi_comm, CaptureError(function));
  if (failure.has_value())
    throw std::runtime_error("Rank " + std::to_string(failure->first) + ": " + failure->second);
}

} // namespace opensn
