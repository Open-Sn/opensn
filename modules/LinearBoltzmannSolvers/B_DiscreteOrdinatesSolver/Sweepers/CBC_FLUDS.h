#ifndef CHITECH_CBC_FLUDS_H
#define CHITECH_CBC_FLUDS_H

#include "mesh/SweepUtilities/FLUDS/FLUDS.h"
#include "CBC_FLUDSCommonData.h"

#include <map>
#include <functional>

namespace chi_math
{
class UnknownManager;
class SpatialDiscretization;
} // namespace chi_math

namespace chi_mesh
{
class Cell;
}

namespace lbs
{

class CBC_FLUDS : public chi_mesh::sweep_management::FLUDS
{
public:
  CBC_FLUDS(size_t num_groups,
            size_t num_angles,
            const CBC_FLUDSCommonData& common_data,
            std::vector<double>& local_psi_data,
            const chi_math::UnknownManager& psi_uk_man,
            const chi_math::SpatialDiscretization& sdm);

  const chi_mesh::sweep_management::FLUDSCommonData& CommonData() const;

  const std::vector<double>& GetLocalUpwindDataBlock() const;

  const double* GetLocalCellUpwindPsi(const std::vector<double>& psi_data_block,
                                      const chi_mesh::Cell& cell);

  const std::vector<double>& GetNonLocalUpwindData(uint64_t cell_global_id,
                                                   unsigned int face_id) const;

  const double* GetNonLocalUpwindPsi(const std::vector<double>& psi_data,
                                     unsigned int face_node_mapped,
                                     unsigned int angle_set_index);

  void ClearLocalAndReceivePsi() override
  {
    deplocs_outgoing_messages_.clear();
  }
  void ClearSendPsi() override {}
  void AllocateInternalLocalPsi(size_t num_grps, size_t num_angles) override {}
  void AllocateOutgoingPsi(size_t num_grps,
                           size_t num_angles,
                           size_t num_loc_sucs) override
  {
  }

  void AllocateDelayedLocalPsi(size_t num_grps, size_t num_angles) override {}
  void AllocatePrelocIOutgoingPsi(size_t num_grps,
                                  size_t num_angles,
                                  size_t num_loc_deps) override
  {
  }
  void AllocateDelayedPrelocIOutgoingPsi(size_t num_grps,
                                         size_t num_angles,
                                         size_t num_loc_deps) override
  {
  }

  std::vector<double>& DelayedLocalPsi() override { return delayed_local_psi_; }
  std::vector<double>& DelayedLocalPsiOld() override
  {
    return delayed_local_psi_old_;
  }

  std::vector<std::vector<double>>& DeplocIOutgoingPsi() override
  {
    return deplocI_outgoing_psi_;
  }

  std::vector<std::vector<double>>& PrelocIOutgoingPsi() override
  {
    return prelocI_outgoing_psi_;
  }

  std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsi() override
  {
    return delayed_prelocI_outgoing_psi_;
  }
  std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsiOld() override
  {
    return delayed_prelocI_outgoing_psi_old_;
  }

  // cell_global_id
  // face_id
  typedef std::pair<uint64_t, unsigned int> CellFaceKey;

  std::map<CellFaceKey, std::vector<double>>& DeplocsOutgoingMessages()
  {
    return deplocs_outgoing_messages_;
  }

private:
  const CBC_FLUDSCommonData& common_data_;
  std::reference_wrapper<std::vector<double>> local_psi_data_;
  const chi_math::UnknownManager& psi_uk_man_;
  const chi_math::SpatialDiscretization& sdm_;

  std::vector<double> delayed_local_psi_;
  std::vector<double> delayed_local_psi_old_;
  std::vector<std::vector<double>> deplocI_outgoing_psi_;
  std::vector<std::vector<double>> prelocI_outgoing_psi_;
  std::vector<std::vector<double>> boundryI_incoming_psi_;

  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_;
  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_old_;

  std::map<CellFaceKey, std::vector<double>> deplocs_outgoing_messages_;
};

} // namespace lbs

#endif // CHITECH_CBC_FLUDS_H
