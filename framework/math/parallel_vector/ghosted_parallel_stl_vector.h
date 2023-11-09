#pragma once

#include "framework/math/parallel_vector/parallel_stl_vector.h"

#include "framework/math/vector_ghost_communicator/vector_ghost_communicator.h"

namespace chi_math
{

class GhostedParallelSTLVector : public ParallelSTLVector
{
public:
  /**
   * Initialize the ghosted parallel vector with the given local and global
   * sizes, with the specified global ghost indices.
   */
  GhostedParallelSTLVector(const uint64_t local_size,
                           const uint64_t global_size,
                           const std::vector<int64_t>& ghost_ids,
                           const MPI_Comm communicator = MPI_COMM_WORLD)
    : ParallelSTLVector(local_size, global_size, communicator),
      ghost_comm_(local_size, global_size, ghost_ids, communicator)
  {
    values_.assign(local_size_ + ghost_comm_.NumGhosts(), 0.0);
  }

  /// Initialize a ghosted parallel vector from a ghost communicator.
  explicit GhostedParallelSTLVector(const VectorGhostCommunicator& ghost_comm)
    : ParallelSTLVector(ghost_comm.LocalSize(), ghost_comm.GlobalSize(), ghost_comm.Communicator()),
      ghost_comm_(ghost_comm)
  {
    values_.assign(local_size_ + ghost_comm_.NumGhosts(), 0.0);
  }

  /**Copy constructor.*/
  GhostedParallelSTLVector(const GhostedParallelSTLVector& other)
    : ParallelSTLVector(other), ghost_comm_(other.ghost_comm_)
  {
  }

  /**Move constructor.*/
  GhostedParallelSTLVector(GhostedParallelSTLVector&& other) noexcept
    : ParallelSTLVector(std::move(other)), ghost_comm_(std::move(other.ghost_comm_))
  {
  }

  std::unique_ptr<ParallelVector> MakeCopy() const override;
  std::unique_ptr<ParallelVector> MakeClone() const override;

  /// Return the number of ghosts associated with the local vector.
  uint64_t NumGhosts() const { return ghost_comm_.NumGhosts(); }

  /// Return the total size of the local vector, including ghosts.
  uint64_t LocalSizeWithGhosts() const { return values_.size(); }

  /// Return the ghost indices associated with the local vector.
  const std::vector<int64_t>& GhostIndices() const { return ghost_comm_.GhostIndices(); }

  /// Map a global ghost id to its respective local id.
  int64_t MapGhostToLocal(const int64_t ghost_id) const
  {
    return ghost_comm_.MapGhostToLocal(ghost_id);
  }

  /// Return a vector containing the locally owned data and ghost data.
  std::vector<double> MakeGhostedLocalVector() const { return values_; }

  /**
   * Return the value of the parallel vector for the specified global index.
   *
   * An error is thrown if the global index does not belong to a locally
   * owned, or ghost entry.
   */
  double GetGlobalValue(int64_t global_id) const;

  /**
   * Communicate the current ghost entries to all other processes to
   * update the locally stored ghost data.
   */
  void CommunicateGhostEntries() override { ghost_comm_.CommunicateGhostEntries(values_); }

private:
  VectorGhostCommunicator ghost_comm_;
};

} // namespace chi_math
