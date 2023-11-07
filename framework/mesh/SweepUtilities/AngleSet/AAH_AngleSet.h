#pragma once

#include "framework/mesh/SweepUtilities/AngleSet/AngleSet.h"

namespace chi_mesh::sweep_management
{

/**Manages the workstages of a single angle set.*/
class AAH_AngleSet : public AngleSet
{
public:
  /**AngleSet constructor.*/
  AAH_AngleSet(size_t id,
               size_t in_numgrps,
               size_t in_ref_subset,
               const SPDS& in_spds,
               std::shared_ptr<FLUDS>& in_fluds,
               std::vector<size_t>& angle_indices,
               std::map<uint64_t, std::shared_ptr<SweepBndry>>& sim_boundaries,
               int sweep_eager_limit,
               const chi::ChiMPICommunicatorSet& in_comm_set);

  void InitializeDelayedUpstreamData() override;

  int GetMaxBufferMessages() const override;

  void SetMaxBufferMessages(int new_max) override;

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk,
                                 const std::vector<size_t>& timing_tags,
                                 ExecutionPermission permission) override;
  AngleSetStatus FlushSendBuffers() override;
  void ResetSweepBuffers() override;
  bool ReceiveDelayedData() override;

  const double* PsiBndry(uint64_t bndry_map,
                         unsigned int angle_num,
                         uint64_t cell_local_id,
                         unsigned int face_num,
                         unsigned int fi,
                         int g,
                         size_t gs_ss_begin,
                         bool surface_source_active) override;
  double* ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                                     unsigned int angle_num,
                                     uint64_t cell_local_id,
                                     unsigned int face_num,
                                     unsigned int fi,
                                     size_t gs_ss_begin) override;

protected:
  chi_mesh::sweep_management::AAH_ASynchronousCommunicator async_comm_;
};

} // namespace chi_mesh::sweep_management
