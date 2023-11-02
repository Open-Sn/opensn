#ifndef CHITECH_CBC_ANGLESET_H
#define CHITECH_CBC_ANGLESET_H

#include "mesh/SweepUtilities/AngleSet/AngleSet.h"
#include "CBC_AsyncComm.h"

namespace lbs
{

class CBC_SPDS;
struct Task;

class CBC_AngleSet : public chi_mesh::sweep_management::AngleSet
{
public:
  CBC_AngleSet(size_t id,
               size_t num_groups,
               const chi_mesh::sweep_management::SPDS& spds,
               std::shared_ptr<chi_mesh::sweep_management::FLUDS>& fluds,
               const std::vector<size_t>& angle_indices,
               std::map<uint64_t, SweepBndryPtr>& sim_boundaries,
               size_t in_ref_subset,
               const chi::ChiMPICommunicatorSet& comm_set);

  chi_mesh::sweep_management::AsynchronousCommunicator*
  GetCommunicator() override;
  void InitializeDelayedUpstreamData() override {}
  int GetMaxBufferMessages() const override { return 0; }
  void SetMaxBufferMessages(int new_max) override {}

  chi_mesh::sweep_management::AngleSetStatus AngleSetAdvance(
    chi_mesh::sweep_management::SweepChunk& sweep_chunk,
    const std::vector<size_t>& timing_tags,
    chi_mesh::sweep_management::ExecutionPermission permission) override;

  chi_mesh::sweep_management::AngleSetStatus FlushSendBuffers() override
  {
    const bool all_messages_sent = async_comm_.SendData();
    return all_messages_sent
             ? chi_mesh::sweep_management::AngleSetStatus::MESSAGES_SENT
             : chi_mesh::sweep_management::AngleSetStatus::MESSAGES_PENDING;
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
  std::vector<chi_mesh::sweep_management::Task> current_task_list_;
  CBC_ASynchronousCommunicator async_comm_;
};

} // namespace lbs

#endif // CHITECH_CBC_ANGLESET_H
