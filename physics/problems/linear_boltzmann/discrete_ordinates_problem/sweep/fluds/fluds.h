// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <vector>
#include <set>
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
      spds_(spds){};

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

  virtual std::vector<double>& DelayedLocalPsi() = 0;
  virtual std::vector<double>& DelayedLocalPsiOld() = 0;

  virtual std::vector<std::vector<double>>& DeplocIOutgoingPsi() = 0;

  virtual std::vector<std::vector<double>>& PrelocIOutgoingPsi() = 0;

  virtual std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsi() = 0;

  virtual std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsiOld() = 0;

  virtual ~FLUDS() = default;

protected:
  const size_t num_groups_;
  const size_t num_angles_;
  const size_t num_groups_and_angles_;
  const SPDS& spds_;
};

} // namespace opensn
