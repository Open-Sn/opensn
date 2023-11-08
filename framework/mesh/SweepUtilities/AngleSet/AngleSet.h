#pragma once

#include "framework/mesh/chi_mesh.h"
#include "framework/mesh/SweepUtilities/Communicators/AAH_AsynComm.h"
#include "framework/mesh/SweepUtilities/SweepBoundary/sweep_boundary.h"
#include "framework/mesh/SweepUtilities/FLUDS/FLUDS.h"

#include "framework/mpi/chi_mpi.h"

typedef chi_mesh::sweep_management::SweepBoundary SweepBndry;

#include <memory>

namespace chi_mesh::sweep_management
{

class AngleSet
{
public:
  typedef std::shared_ptr<SweepBndry> SweepBndryPtr;

  /**AngleSet constructor.*/
  AngleSet(size_t id,
           size_t num_groups,
           const SPDS& spds,
           std::shared_ptr<FLUDS>& fluds,
           const std::vector<size_t>& angle_indices,
           std::map<uint64_t, SweepBndryPtr>& sim_boundaries,
           size_t in_ref_subset);

  /**Returns the angleset's unique id.*/
  size_t GetID() const;
  /**Returns a reference to the associated spds.*/
  const SPDS& GetSPDS() const;
  /**Returns a reference to the associated fluds_.*/
  FLUDS& GetFLUDS();
  /**Return the reference group subset number.*/
  size_t GetRefGroupSubset() const;
  /**Returns the angle indices associated with this angleset.*/
  const std::vector<size_t>& GetAngleIndices() const;
  /**Returns the angle indices associated with this angleset.*/
  std::map<uint64_t, SweepBndryPtr>& GetBoundaries();

  size_t GetNumGroups() const;
  size_t GetNumAngles() const;

  // Virtual methods
  virtual AsynchronousCommunicator* GetCommunicator();
  /**Initializes delayed upstream data. This method gets called
   * when a sweep scheduler is constructed.*/
  virtual void InitializeDelayedUpstreamData() = 0;

  /**Returns the maximum buffer size from the sweepbuffer.*/
  virtual int GetMaxBufferMessages() const = 0;

  /**Sets the maximum buffer size for the sweepbuffer.*/
  virtual void SetMaxBufferMessages(int new_max) = 0;

  /**This function advances the work stages of an angleset.*/
  virtual AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk,
                                         const std::vector<size_t>& timing_tags,
                                         ExecutionPermission permission) = 0;
  virtual AngleSetStatus FlushSendBuffers() = 0;
  /**Resets the sweep buffer.*/
  virtual void ResetSweepBuffers() = 0;
  /**Instructs the sweep buffer to receive delayed data.*/
  virtual bool ReceiveDelayedData() = 0;

  /**Returns a pointer to a boundary flux data.*/
  virtual const double* PsiBndry(uint64_t bndry_map,
                                 unsigned int angle_num,
                                 uint64_t cell_local_id,
                                 unsigned int face_num,
                                 unsigned int fi,
                                 int g,
                                 size_t gs_ss_begin,
                                 bool surface_source_active) = 0;
  /**Returns a pointer to outbound boundary flux data.*/
  virtual double* ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                                             unsigned int angle_num,
                                             uint64_t cell_local_id,
                                             unsigned int face_num,
                                             unsigned int fi,
                                             size_t gs_ss_begin) = 0;

  virtual ~AngleSet() = default;

protected:
  const size_t id_;
  const size_t num_grps;
  const SPDS& spds_;
  std::shared_ptr<FLUDS> fluds_;
  const std::vector<size_t> angles_;
  std::map<uint64_t, SweepBndryPtr>& ref_boundaries_;
  const size_t ref_group_subset_;

  bool executed_ = false;
};

} // namespace chi_mesh::sweep_management
