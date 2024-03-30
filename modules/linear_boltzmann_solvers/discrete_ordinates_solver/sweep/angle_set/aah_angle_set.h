#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/communicators/aah_async_comm.h"

namespace opensn
{
namespace lbs
{

/**Manages the workstages of a single angle set.*/
class AAH_AngleSet : public AngleSet
{
protected:
  AAH_ASynchronousCommunicator async_comm_;

public:
  AAH_AngleSet(size_t id,
               size_t num_groups,
               size_t group_subset,
               const SPDS& spds,
               std::shared_ptr<FLUDS>& fluds,
               std::vector<size_t>& angle_indices,
               std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
               int maximum_message_size,
               const MPICommunicatorSet& in_comm_set);

  void InitializeDelayedUpstreamData() override;

  int GetMaxBufferMessages() const override;

  void SetMaxBufferMessages(int new_max) override;

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk,
                                 AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override;

  void ResetSweepBuffers() override;

  bool ReceiveDelayedData() override;

  const double* PsiBoundary(uint64_t boundary_id,
                            unsigned int angle_num,
                            uint64_t cell_local_id,
                            unsigned int face_num,
                            unsigned int fi,
                            int g,
                            size_t gs_ss_begin,
                            bool surface_source_active) override;

  double* PsiReflected(uint64_t boundary_id,
                       unsigned int angle_num,
                       uint64_t cell_local_id,
                       unsigned int face_num,
                       unsigned int fi,
                       size_t gs_ss_begin) override;
};

} // namespace lbs
} // namespace opensn
