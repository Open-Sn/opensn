#pragma once

#include "opensn/framework/mesh/chi_mesh.h"
#include "opensn/framework/mesh/SweepUtilities/Communicators/AAH_AsynComm.h"
#include "opensn/framework/mesh/SweepUtilities/SweepBoundary/sweep_boundary.h"
#include "opensn/framework/mesh/SweepUtilities/FLUDS/FLUDS.h"

#include "opensn/framework/mpi/chi_mpi.h"

typedef chi_mesh::sweep_management::SweepBoundary SweepBndry;

#include <memory>

namespace chi_mesh::sweep_management
{

class AngleSet
{
public:
  typedef std::shared_ptr<SweepBndry> SweepBndryPtr;

  AngleSet(size_t id,
           size_t num_groups,
           const SPDS& spds,
           std::shared_ptr<FLUDS>& fluds,
           const std::vector<size_t>& angle_indices,
           std::map<uint64_t, SweepBndryPtr>& sim_boundaries,
           size_t in_ref_subset);

  size_t GetID() const;
  const SPDS& GetSPDS() const;
  FLUDS& GetFLUDS();
  size_t GetRefGroupSubset() const;
  const std::vector<size_t>& GetAngleIndices() const;
  std::map<uint64_t, SweepBndryPtr>& GetBoundaries();

  size_t GetNumGroups() const;
  size_t GetNumAngles() const;

  // Virtual methods
  virtual AsynchronousCommunicator* GetCommunicator();
  virtual void InitializeDelayedUpstreamData() = 0;

  virtual int GetMaxBufferMessages() const = 0;

  virtual void SetMaxBufferMessages(int new_max) = 0;

  virtual AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk,
                                         const std::vector<size_t>& timing_tags,
                                         ExecutionPermission permission) = 0;
  virtual AngleSetStatus FlushSendBuffers() = 0;
  virtual void ResetSweepBuffers() = 0;
  virtual bool ReceiveDelayedData() = 0;

  virtual const double* PsiBndry(uint64_t bndry_map,
                                 unsigned int angle_num,
                                 uint64_t cell_local_id,
                                 unsigned int face_num,
                                 unsigned int fi,
                                 int g,
                                 size_t gs_ss_begin,
                                 bool surface_source_active) = 0;
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
