#pragma once

#include "framework/mesh/sweep_utilities/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_async_comm.h"

namespace opensn
{
struct Task;

namespace lbs
{
class CBC_SPDS;

class CBC_AngleSet : public AngleSet
{
public:
  CBC_AngleSet(size_t id,
               size_t num_groups,
               const SPDS& spds,
               std::shared_ptr<FLUDS>& fluds,
               const std::vector<size_t>& angle_indices,
               std::map<uint64_t, SweepBndryPtr>& sim_boundaries,
               size_t in_ref_subset,
               const ChiMPICommunicatorSet& comm_set);

  AsynchronousCommunicator* GetCommunicator() override;
  void InitializeDelayedUpstreamData() override {}
  int GetMaxBufferMessages() const override { return 0; }
  void SetMaxBufferMessages(int new_max) override {}

  AngleSetStatus AngleSetAdvance(opensn::SweepChunk& sweep_chunk,
                                 const std::vector<size_t>& timing_tags,
                                 ExecutionPermission permission) override;

  AngleSetStatus FlushSendBuffers() override
  {
    const bool all_messages_sent = async_comm_.SendData();
    return all_messages_sent ? AngleSetStatus::MESSAGES_SENT : AngleSetStatus::MESSAGES_PENDING;
  }
  void ResetSweepBuffers() override;
  bool ReceiveDelayedData() override { return true; }
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
  const CBC_SPDS& cbc_spds_;
  std::vector<Task> current_task_list_;
  CBC_ASynchronousCommunicator async_comm_;
};

} // namespace lbs
} // namespace opensn
