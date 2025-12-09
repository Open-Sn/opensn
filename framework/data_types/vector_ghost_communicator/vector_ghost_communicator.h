// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpicpp-lite/mpicpp-lite.h"
#include <vector>
#include <cstdint>
#include <map>
#include <optional>

namespace mpi = mpicpp_lite;

namespace opensn
{

/// Vector with allocation space for ghosts.
class VectorGhostCommunicator
{
public:
  VectorGhostCommunicator(uint64_t local_size,
                          uint64_t global_size,
                          const std::vector<uint64_t>& ghost_ids,
                          const mpi::Communicator& communicator);

  /// Copy constructor.
  VectorGhostCommunicator(const VectorGhostCommunicator& other);

  /// Move constructor.
  VectorGhostCommunicator(VectorGhostCommunicator&& other) = default;

  uint64_t GetLocalSize() const { return local_size_; }
  uint64_t GetGlobalSize() const { return global_size_; }

  uint64_t GetNumGhosts() const { return ghost_ids_.size(); }
  const std::vector<uint64_t>& GetGhostIndices() const { return ghost_ids_; }

  const mpi::Communicator& GetCommunicator() const { return comm_; }

  std::optional<uint64_t> MapGhostToLocal(uint64_t ghost_id) const;

  void CommunicateGhostEntries(std::vector<double>& ghosted_vector) const;

  std::vector<double> MakeGhostedVector() const;
  std::vector<double> MakeGhostedVector(const std::vector<double>& local_vector) const;

protected:
  const uint64_t local_size_;
  const uint64_t global_size_;
  const std::vector<uint64_t> ghost_ids_;
  const mpi::Communicator& comm_;

  const int location_id_;
  const int process_count_;
  const std::vector<uint64_t> extents_;

  struct CachedParallelData
  {
    std::vector<int> sendcounts;
    std::vector<int> senddispls;
    std::vector<int> recvcounts;
    std::vector<int> recvdispls;

    std::vector<uint64_t> local_ids_to_send;
    std::map<uint64_t, size_t> ghost_to_recv_map;
  };

  const CachedParallelData cached_parallel_data_;

private:
  int FindOwnerPID(uint64_t global_id) const;
  CachedParallelData MakeCachedParallelData();
};

} // namespace opensn
