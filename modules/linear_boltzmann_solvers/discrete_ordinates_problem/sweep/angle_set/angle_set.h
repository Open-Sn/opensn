// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "framework/mesh/mesh.h"
#include "framework/logging/log.h"
#include <memory>

namespace opensn
{

class SweepChunk;

class AngleSet
{
public:
  AngleSet(size_t id,
           const LBSGroupset& groupset,
           const SPDS& spds,
           std::shared_ptr<FLUDS>& fluds,
           const std::vector<size_t>& angle_indices,
           std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries);

  /// Returns the angleset's unique id.
  size_t GetID() const { return id_; }

  /// Returns the associated groupset id.
  int GetGroupsetID() const { return groupset_id_; }

  /// Returns a reference to the associated spds.
  const SPDS& GetSPDS() const { return spds_; }

  /// Returns a reference to the associated fluds_.
  FLUDS& GetFLUDS() { return *fluds_; }

  /// Returns the angle indices associated with this angleset.
  std::vector<std::uint32_t>& GetAngleIndices() { return angles_; }

  /// Returns the angle indices associated with this angleset.
  const std::vector<std::uint32_t>& GetAngleIndices() const { return angles_; }

  /// Returns the angle indices associated with this angleset.
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetBoundaries() { return boundaries_; }

  unsigned int GetNumGroups() const { return num_groups_; }

  size_t GetNumAngles() const { return angles_.size(); }

  bool HasAngleIndex(std::uint32_t angle_index) const;

  /// Update following angle sets.
  void UpdateSweepDependencies(std::set<AngleSet*>& following_angle_sets);

  virtual AsynchronousCommunicator* GetCommunicator()
  {
    OpenSnLogicalError("Method not implemented");
  }

  /**
   * Initializes delayed upstream data. This method gets called when a sweep scheduler is
   * constructed.
   */
  virtual void InitializeDelayedUpstreamData() = 0;

  /// Returns the maximum buffer size from the sweepbuffer.
  virtual int GetMaxBufferMessages() const = 0;

  /// Sets the maximum buffer size for the sweepbuffer.
  virtual void SetMaxBufferMessages(int new_max) = 0;

  /// This function advances the work stages of an angleset.
  virtual AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) = 0;

  virtual AngleSetStatus FlushSendBuffers() = 0;

  /// Resets the sweep buffer.
  virtual void ResetSweepBuffers() = 0;

  /// Resets dependency counter.
  void ResetDependencyCounter() { dependency_counter_ = num_dependencies_; }

  /// Checks if dependency is resolved.
  bool IsDependencyResolved() const { return dependency_counter_ == 0; }

  void DecrementCounter() { --dependency_counter_; }

  /// Instructs the sweep buffer to receive delayed data.
  virtual bool ReceiveDelayedData() = 0;

  /// Returns a pointer to a boundary flux data.
  const double* PsiBoundary(uint64_t boundary_id,
                            unsigned int angle_num,
                            uint64_t cell_local_id,
                            unsigned int face_num,
                            unsigned int fi,
                            unsigned int g,
                            bool surface_source_active)
  {
    const auto& boundary = boundaries_[boundary_id];
    if (boundary->IsReflecting() || surface_source_active)
      return boundary->PsiIncoming(cell_local_id, face_num, fi, angle_num, groupset_id_, g);
    return boundary->ZeroFlux(groupset_id_, g);
  }

  /// Returns a pointer to outbound reflected flux data.
  double* PsiReflected(uint64_t boundary_id,
                       unsigned int angle_num,
                       uint64_t cell_local_id,
                       unsigned int face_num,
                       unsigned int fi)
  {
    return boundaries_[boundary_id]->PsiOutgoing(
      cell_local_id, face_num, fi, angle_num, groupset_id_);
  }

  /// Update the angle index data on the device to match the host.
  virtual void SyncDeviceAngleIndices() {}

  virtual ~AngleSet() = default;

protected:
  const size_t id_;
  const int groupset_id_;
  const unsigned int num_groups_;
  const SPDS& spds_;
  std::shared_ptr<FLUDS> fluds_;
  std::vector<std::uint32_t> angles_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries_;
  bool executed_ = false;

  /// Number of anglesets the current angle set depends on.
  std::size_t num_dependencies_ = 0;
  /// Counter for un-resolved dependencies.
  std::size_t dependency_counter_ = 0;
  /**
   * List of angle sets waiting after this angle set.
   * After this angle set completes its sweep chunk, it decrements the counter of the angle sets
   * waiting after it to allow them to proceed.
   */
  std::vector<AngleSet*> following_angle_sets_;
};

} // namespace opensn
