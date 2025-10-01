// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <set>
#include <span>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace opensn
{

class GridFaceHistogram;
class SPDS;

class FLUDS
{
public:
  FLUDS(size_t num_groups, size_t num_angles, const SPDS& spds)
    : num_groups_(num_groups),
      num_angles_(num_angles),
      num_groups_and_angles_(num_groups_ * num_angles_),
      spds_(spds) {};

  const SPDS& GetSPDS() const { return spds_; }

  virtual void ClearLocalAndReceivePsi() {}
  virtual void ClearSendPsi() {}
  virtual void AllocateInternalLocalPsi(size_t num_grps, size_t num_angles) {}
  virtual void AllocateOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_sucs) {}

  virtual void AllocateDelayedLocalPsi(size_t num_grps, size_t num_angles) {}
  virtual void AllocatePrelocIOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_deps)
  {
  }
  virtual void
  AllocateDelayedPrelocIOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_deps)
  {
  }

  std::span<double>& DelayedLocalPsi() { return delayed_local_psi_view_; }
  std::span<double>& DelayedLocalPsiOld() { return delayed_local_psi_old_view_; }
  virtual void SetDelayedLocalPsiOldToNew() {}
  virtual void SetDelayedLocalPsiNewToOld() {}

  std::vector<std::span<double>>& DeplocIOutgoingPsi() { return deplocI_outgoing_psi_view_; }

  std::vector<std::span<double>>& PrelocIOutgoingPsi() { return prelocI_outgoing_psi_view_; }

  std::vector<std::span<double>>& DelayedPrelocIOutgoingPsi()
  {
    return delayed_prelocI_outgoing_psi_view_;
  }
  std::vector<std::span<double>>& DelayedPrelocIOutgoingPsiOld()
  {
    return delayed_prelocI_outgoing_psi_old_view_;
  }
  virtual void SetDelayedOutgoingPsiOldToNew() {}
  virtual void SetDelayedOutgoingPsiNewToOld() {}

  virtual ~FLUDS() = default;

protected:
  const size_t num_groups_;
  const size_t num_angles_;
  const size_t num_groups_and_angles_;
  const SPDS& spds_;

  std::span<double> delayed_local_psi_view_;
  std::span<double> delayed_local_psi_old_view_;
  std::vector<std::span<double>> deplocI_outgoing_psi_view_;
  std::vector<std::span<double>> prelocI_outgoing_psi_view_;
  std::vector<std::span<double>> delayed_prelocI_outgoing_psi_view_;
  std::vector<std::span<double>> delayed_prelocI_outgoing_psi_old_view_;
};

} // namespace opensn
