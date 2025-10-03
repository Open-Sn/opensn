// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/parallel_vector/parallel_stl_vector.h"
#include "framework/data_types/vector_ghost_communicator/vector_ghost_communicator.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <stdexcept>

namespace mpi = mpicpp_lite;

namespace opensn
{

class GhostedParallelSTLVector : public ParallelSTLVector
{
public:
  /**
   * Initialize the ghosted parallel vector with the given local and global sizes, with the
   * specified global ghost indices.
   */
  GhostedParallelSTLVector(uint64_t local_size,
                           uint64_t global_size,
                           const std::vector<int64_t>& ghost_ids,
                           const mpi::Communicator& communicator = MPI_COMM_WORLD)
    : ParallelSTLVector(local_size, global_size, communicator),
      ghost_comm_(local_size, global_size, ghost_ids, communicator)
  {
    values_.assign(local_size_ + ghost_comm_.GetNumGhosts(), 0.0);
  }

  /// Initialize a ghosted parallel vector from a ghost communicator.
  explicit GhostedParallelSTLVector(const VectorGhostCommunicator& ghost_comm)
    : ParallelSTLVector(
        ghost_comm.GetLocalSize(), ghost_comm.GetGlobalSize(), ghost_comm.GetCommunicator()),
      ghost_comm_(ghost_comm)
  {
    values_.assign(local_size_ + ghost_comm_.GetNumGhosts(), 0.0);
  }

  /// Copy constructor.
  GhostedParallelSTLVector(const GhostedParallelSTLVector& other) = default;

  /// Move constructor.
  GhostedParallelSTLVector(GhostedParallelSTLVector&& other) noexcept = default;

  std::unique_ptr<ParallelVector> MakeCopy() const override;
  std::unique_ptr<ParallelVector> MakeClone() const override;

  /// Return the number of ghosts associated with the local vector.
  uint64_t GetNumGhosts() const { return ghost_comm_.GetNumGhosts(); }

  /// Return the total size of the local vector, including ghosts.
  uint64_t GetLocalSizeWithGhosts() const { return values_.size(); }

  /// Return the ghost indices associated with the local vector.
  const std::vector<int64_t>& GetGhostIndices() const { return ghost_comm_.GetGhostIndices(); }

  /// Map a global ghost id to its respective local id.
  int64_t MapGhostToLocal(const int64_t ghost_id) const
  {
    auto local_id = ghost_comm_.MapGhostToLocal(ghost_id);
    if (local_id)
      return local_id.value();
    else
      throw std::runtime_error("Mapping from ghost ID to local ID failed");
  }

  /// Return a vector containing the locally owned data and ghost data.
  std::vector<double> MakeGhostedLocalVector() const { return values_; }

  /**
   * Return the value of the parallel vector for the specified global index.
   *
   * An error is thrown if the global index does not belong to a locally owned, or ghost entry.
   */
  double GetGlobalValue(uint64_t global_id) const;

  /**
   * Communicate the current ghost entries to all other processes to update the locally stored ghost
   * data.
   */
  void CommunicateGhostEntries() override { ghost_comm_.CommunicateGhostEntries(values_); }

private:
  VectorGhostCommunicator ghost_comm_;
};

} // namespace opensn
