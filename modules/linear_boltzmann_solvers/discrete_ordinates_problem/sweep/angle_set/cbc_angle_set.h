// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"

namespace opensn
{

class CBC_FLUDS;
class CBC_SPDS;

/**
 * Host-side CBC angle set.
 *
 * Owns the local CBC task state for one angle set and advances the host CBC sweep
 * by combining local task execution with non-local message progress.
 */
class CBC_AngleSet : public AngleSet
{
public:
  /**
   * Construct the CBC angle set.
   *
   * \param id Angle-set ID.
   * \param num_groups Number of groups in the angle set.
   * \param spds Sweep plane data structure for this angle set.
   * \param fluds CBC FLUDS instance for this angle set.
   * \param angle_indices Global angle indices represented by this angle set.
   * \param boundaries Sweep-boundary table indexed by boundary ID.
   * \param comm_set MPI communicator set used for receives.
   */
  CBC_AngleSet(size_t id,
               unsigned int num_groups,
               const SPDS& spds,
               std::shared_ptr<FLUDS>& fluds,
               const std::vector<size_t>& angle_indices,
               std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
               const MPICommunicatorSet& comm_set);

  /// Return the delayed-data communicator for this angle set.
  AsynchronousCommunicator* GetCommunicator() override;

  /// Initialize delayed upstream data before the sweep starts.
  void InitializeDelayedUpstreamData() override {}

  /// Return the buffered-message limit used by the scheduler.
  int GetMaxBufferMessages() const override { return 0; }

  /// Set the buffered-message limit used by the scheduler.
  void SetMaxBufferMessages(int max_buffer_messages) override {}

  /// Advance the host CBC angle set by one scheduler step.
  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override
  {
    const bool all_messages_sent = async_comm_.SendData();
    return all_messages_sent ? AngleSetStatus::MESSAGES_SENT : AngleSetStatus::MESSAGES_PENDING;
  }

  void ResetSweepBuffers() override;

  /// Report whether delayed upstream data has been received.
  bool ReceiveDelayedData() override { return true; }

  /**
   * Return the incoming boundary angular flux for one boundary face node.
   *
   * \param boundary_id Boundary ID.
   * \param angle_num Angle index within the angle set.
   * \param cell_local_id Local cell ID.
   * \param face_num Face ID on the cell.
   * \param fi Face-node index.
   * \param g Group index.
   * \param surface_source_active Flag if surface source is active.
   * \return Pointer to the requested incoming boundary value.
   */
  const double* PsiBoundary(uint64_t boundary_id,
                            unsigned int angle_num,
                            uint64_t cell_local_id,
                            unsigned int face_num,
                            unsigned int fi,
                            unsigned int g,
                            bool surface_source_active) override;

  /**
   * Return the outgoing reflecting-boundary storage for one face node.
   *
   * \param boundary_id Boundary ID.
   * \param angle_num Angle index within the angle set.
   * \param cell_local_id Local cell ID.
   * \param face_num Face ID on the cell.
   * \param fi Face-node index.
   * \return Pointer to the reflected outgoing storage for the node.
   */
  double* PsiReflected(uint64_t boundary_id,
                       unsigned int angle_num,
                       uint64_t cell_local_id,
                       unsigned int face_num,
                       unsigned int fi) override;

protected:
  /// Reset the mutable local-task state before a new sweep.
  void ResetTaskState();

  /// CBC sweep plane data structure for this angle set.
  const CBC_SPDS& cbc_spds_;
  /// Initial predecessor counts per local CBC task.
  std::vector<unsigned int> initial_dependencies_;
  /// Local tasks that are ready at the start of a sweep.
  std::vector<std::uint32_t> initial_ready_tasks_;
  /// Mutable predecessor counts for the current sweep.
  std::vector<unsigned int> remaining_dependencies_;
  /// Local tasks ready to execute.
  std::vector<std::uint32_t> ready_tasks_;
  /// Number of completed local tasks.
  size_t num_completed_tasks_ = 0;
  /// Asynchronous communicator for this angle set.
  CBC_AsynchronousCommunicator async_comm_;
  /// CBC FLUDS instance for this angle set.
  CBC_FLUDS& cbc_fluds_;
};

} // namespace opensn
