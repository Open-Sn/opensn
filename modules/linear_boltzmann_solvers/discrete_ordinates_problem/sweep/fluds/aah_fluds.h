// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"

namespace opensn
{

/// Implementation of the Adams-Adams-Hawkins Flux data structure.
class AAH_FLUDS : public FLUDS
{
public:
  /**
   * This constructor initializes an auxiliary FLUDS based on a primary FLUDS. The restriction here
   * is that the auxiliary FLUDS has the exact same sweep ordering as the primary FLUDS.
   */
  AAH_FLUDS(size_t num_groups, size_t num_angles, const AAH_FLUDSCommonData& common_data);

private:
  const AAH_FLUDSCommonData& common_data_;

  // local_psi_n_block_stride[fc]. Given face category fc, the value is
  // total number of faces that store information in this category's buffer
  // per angle
  std::vector<size_t> local_psi_Gn_block_strideG_; // Custom G

  size_t delayed_local_psi_Gn_block_strideG_; // Custom G

  std::vector<std::vector<double>> local_psi_;
  std::vector<double> delayed_local_psi_;
  std::vector<double> delayed_local_psi_old_;
  std::vector<std::vector<double>> deplocI_outgoing_psi_;
  std::vector<std::vector<double>> prelocI_outgoing_psi_;
  std::vector<std::vector<double>> boundryI_incoming_psi_;

  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_;
  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_old_;

public:
  /**
   * Given a sweep ordering index, the outgoing face counter, the outgoing face dof, this function
   * computes the location of this position's upwind psi in the local upwind psi vector and returns
   * a reference to it.
   */
  double* OutgoingPsi(std::size_t cell_so_index,
                      int outb_face_counter,
                      std::size_t face_dof,
                      std::size_t n);

  /**
   * Given a sweep ordering index, the incoming face counter, the incoming face dof, this function
   * computes the location where to store this position's outgoing psi and returns a reference to
   * it.
   */
  double* UpwindPsi(std::size_t cell_so_index,
                    int inc_face_counter,
                    std::size_t face_dof,
                    std::size_t g,
                    std::size_t n);

  /// Given a outbound face counter this method returns a pointer to the location
  double* NLOutgoingPsi(int outb_face_count, std::size_t face_dof, std::size_t n);

  /**
   * Given a sweep ordering index, the incoming face counter, the incoming face dof, this function
   * computes the location where to obtain the position's upwind psi.
   */
  double*
  NLUpwindPsi(int nonl_inc_face_counter, std::size_t face_dof, std::size_t g, std::size_t n);

  size_t GetPrelocIFaceDOFCount(int prelocI) const;
  size_t GetDelayedPrelocIFaceDOFCount(int prelocI) const;
  size_t GetDeplocIFaceDOFCount(int deplocI) const;

  void ClearLocalAndReceivePsi() override;
  void ClearSendPsi() override;
  void AllocateInternalLocalPsi() override;
  void AllocateOutgoingPsi() override;

  void AllocateDelayedLocalPsi() override;
  void AllocatePrelocIOutgoingPsi() override;
  void AllocateDelayedPrelocIOutgoingPsi() override;

  void SetDelayedOutgoingPsiOldToNew() override;
  void SetDelayedOutgoingPsiNewToOld() override;

  void SetDelayedLocalPsiOldToNew() override;
  void SetDelayedLocalPsiNewToOld() override;
};

} // namespace opensn
