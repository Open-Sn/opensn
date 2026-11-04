// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/runtime.h"
#include <cstdint>
#include <exception>
#include <optional>
#include <stdexcept>
#include <string>
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
 * Find the lowest rank of `comm` whose `local_error` is not empty. Returns that rank and its error
 * on every rank, or nullopt on every rank if no rank failed. Collective.
 */
inline std::optional<std::pair<int, std::string>>
FindFirstRankError(const mpi::Communicator& comm, const std::string& local_error)
{
  const int comm_size = comm.size();
  int failure_rank = local_error.empty() ? comm_size : comm.rank();
  comm.all_reduce(failure_rank, mpi::op::min<int>());
  if (failure_rank == comm_size)
    return std::nullopt;

  std::string message = local_error;
  std::uint64_t length = message.size();
  comm.broadcast(length, failure_rank);
  message.resize(length);
  if (length > 0)
    comm.broadcast(message.data(), static_cast<int>(length), failure_rank);
  return std::make_pair(failure_rank, std::move(message));
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
                     const std::string& operation)
{
  const auto failure = FindFirstRankError(comm, local_error);
  if (failure.has_value())
    throw Exception(operation + " failed on rank " + std::to_string(failure->first) + ": " +
                    failure->second);
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
  std::string error;
  try
  {
    function();
  }
  catch (const std::exception& exception)
  {
    error = DescribeError(exception);
  }
  catch (...)
  {
    error = non_standard_exception_message;
  }
  const auto failure = FindFirstRankError(mpi_comm, error);
  if (failure.has_value())
    throw std::runtime_error("Rank " + std::to_string(failure->first) + ": " + failure->second);
}

/**
 * Run `function`, which makes no MPI calls, on rank `root` of `comm` only, and throw `Exception` on
 * every rank if it throws, as ThrowIfAnyRankFailed does. Collective.
 */
template <typename Exception = std::runtime_error, typename Function>
void
RunOnRootWithCollectiveErrors(const mpi::Communicator& comm,
                              const int root,
                              const Function& function,
                              const std::string& operation)
{
  std::string error;
  if (comm.rank() == root)
  {
    try
    {
      function();
    }
    catch (const std::exception& exception)
    {
      error = DescribeError(exception);
    }
    catch (...)
    {
      error = non_standard_exception_message;
    }
  }
  ThrowIfAnyRankFailed<Exception>(comm, error, operation);
}

} // namespace opensn
