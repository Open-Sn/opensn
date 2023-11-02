#ifndef AUX_FLUDS_h
#define AUX_FLUDS_h

#include "FLUDS.h"
#include "AAH_FLUDSCommonData.h"

namespace chi_mesh::sweep_management
{

// ########################################################################
/**Implementation of the Adams-Adams-Hawkins Flux data structure.*/
class AAH_FLUDS : public FLUDS
{
public:
  AAH_FLUDS(size_t num_groups,
            size_t num_angles,
            const AAH_FLUDSCommonData& common_data);

private:
  const AAH_FLUDSCommonData& common_data_;

  // local_psi_n_block_stride[fc]. Given face category fc, the value is
  // total number of faces that store information in this category's buffer
  // per angle
  std::vector<size_t> local_psi_Gn_block_strideG; // Custom G

  size_t delayed_local_psi_Gn_block_strideG; // Custom G

  std::vector<std::vector<double>> local_psi_;
  std::vector<double> delayed_local_psi_;
  std::vector<double> delayed_local_psi_old_;
  std::vector<std::vector<double>> deplocI_outgoing_psi_;
  std::vector<std::vector<double>> prelocI_outgoing_psi_;
  std::vector<std::vector<double>> boundryI_incoming_psi_;

  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_;
  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_old_;

public:
  double* OutgoingPsi(int cell_so_index,
                      int outb_face_counter,
                      int face_dof,
                      int n);
  double* UpwindPsi(int cell_so_index,
                    int inc_face_counter,
                    int face_dof,
                    int g,
                    int n);

  double* NLOutgoingPsi(int outb_face_count, int face_dof, int n);

  double*
  NLUpwindPsi(int nonl_inc_face_counter, int face_dof, int g, int n);

  size_t GetPrelocIFaceDOFCount(int prelocI) const;
  size_t GetDelayedPrelocIFaceDOFCount(int prelocI) const;
  size_t GetDeplocIFaceDOFCount(int deplocI) const;

  void ClearLocalAndReceivePsi() override;
  void ClearSendPsi() override;
  void AllocateInternalLocalPsi(size_t num_grps, size_t num_angles) override;
  void AllocateOutgoingPsi(size_t num_grps,
                           size_t num_angles,
                           size_t num_loc_sucs) override;

  void AllocateDelayedLocalPsi(size_t num_grps, size_t num_angles) override;
  void AllocatePrelocIOutgoingPsi(size_t num_grps,
                                  size_t num_angles,
                                  size_t num_loc_deps) override;
  void AllocateDelayedPrelocIOutgoingPsi(size_t num_grps,
                                         size_t num_angles,
                                         size_t num_loc_deps) override;

  std::vector<double>& DelayedLocalPsi() override;
  std::vector<double>& DelayedLocalPsiOld() override;

  std::vector<std::vector<double>>& DeplocIOutgoingPsi() override;

  std::vector<std::vector<double>>& PrelocIOutgoingPsi() override;

  std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsi() override;
  std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsiOld() override;
};

} // namespace chi_mesh::sweep_management

#endif