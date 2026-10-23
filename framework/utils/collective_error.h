// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/runtime.h"
#include <exception>
#include <stdexcept>
#include <string>

namespace opensn
{

/**
 * Runs rank-local work that makes no MPI calls and propagates a failure to every rank. If the work
 * throws on any rank, every rank throws std::runtime_error with the message from the lowest failing
 * rank, so no rank is left waiting in a later collective. Collective.
 */
template <typename Function>
void
RunWithCollectiveErrors(const Function& function)
{
  std::string message;
  bool failed = false;
  try
  {
    function();
  }
  catch (const std::exception& error)
  {
    failed = true;
    message = error.what();
  }

  const int comm_size = mpi_comm.size();
  const int local_failure_rank = failed ? mpi_comm.rank() : comm_size;
  int failure_rank = comm_size;
  mpi_comm.all_reduce(local_failure_rank, failure_rank, mpi::op::min<int>());
  if (failure_rank == comm_size)
    return;

  int length = static_cast<int>(message.size());
  mpi_comm.broadcast(length, failure_rank);
  message.resize(length);
  mpi_comm.broadcast(message.data(), length, failure_rank);
  throw std::runtime_error("Rank " + std::to_string(failure_rank) + ": " + message);
}

} // namespace opensn
