#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/sweep.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/fluds.h"
#include "framework/mesh/mesh.h"
#include "framework/logging/log.h"
#include <memory>

namespace opensn
{
namespace lbs
{

class SweepChunk;

class AngleSet
{
protected:
  const size_t id_;
  const size_t num_groups_;
  const SPDS& spds_;
  std::shared_ptr<FLUDS> fluds_;
  const std::vector<size_t> angles_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries_;
  const size_t group_subset_;
  bool executed_ = false;

public:
  AngleSet(size_t id,
           size_t num_groups,
           const SPDS& spds,
           std::shared_ptr<FLUDS>& fluds,
           const std::vector<size_t>& angle_indices,
           std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
           const size_t group_subset)
    : id_(id),
      num_groups_(num_groups),
      spds_(spds),
      fluds_(fluds),
      angles_(angle_indices),
      boundaries_(boundaries),
      group_subset_(group_subset)
  {
  }

  /**Returns the angleset's unique id.*/
  size_t GetID() const { return id_; }

  /**Returns a reference to the associated spds.*/
  const SPDS& GetSPDS() const { return spds_; }

  /**Returns a reference to the associated fluds_.*/
  FLUDS& GetFLUDS() { return *fluds_; }

  /**Return the group subset number.*/
  size_t GetGroupSubset() const { return group_subset_; }

  /**Returns the angle indices associated with this angleset.*/
  const std::vector<size_t>& GetAngleIndices() const { return angles_; }

  /**Returns the angle indices associated with this angleset.*/
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetBoundaries() { return boundaries_; }

  size_t GetNumGroups() const { return num_groups_; }

  size_t GetNumAngles() const { return angles_.size(); }

  virtual AsynchronousCommunicator* GetCommunicator()
  {
    OpenSnLogicalError("Method not implemented");
  }

  /**Initializes delayed upstream data. This method gets called
   * when a sweep scheduler is constructed.*/
  virtual void InitializeDelayedUpstreamData() = 0;

  /**Returns the maximum buffer size from the sweepbuffer.*/
  virtual int GetMaxBufferMessages() const = 0;

  /**Sets the maximum buffer size for the sweepbuffer.*/
  virtual void SetMaxBufferMessages(int new_max) = 0;

  /**This function advances the work stages of an angleset.*/
  virtual AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) = 0;

  virtual AngleSetStatus FlushSendBuffers() = 0;

  /**Resets the sweep buffer.*/
  virtual void ResetSweepBuffers() = 0;

  /**Instructs the sweep buffer to receive delayed data.*/
  virtual bool ReceiveDelayedData() = 0;

  /**Returns a pointer to a boundary flux data.*/
  virtual const double* PsiBoundary(uint64_t boundary_id,
                                    unsigned int angle_num,
                                    uint64_t cell_local_id,
                                    unsigned int face_num,
                                    unsigned int fi,
                                    int g,
                                    size_t gs_ss_begin,
                                    bool surface_source_active) = 0;

  /**Returns a pointer to outbound reflected flux data.*/
  virtual double* PsiReflected(uint64_t boundary_id,
                               unsigned int angle_num,
                               uint64_t cell_local_id,
                               unsigned int face_num,
                               unsigned int fi,
                               size_t gs_ss_begin) = 0;

  virtual ~AngleSet() = default;
};

} // namespace lbs
} // namespace opensn
