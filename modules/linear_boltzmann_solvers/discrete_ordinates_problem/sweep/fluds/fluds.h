// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include "framework/utils/error.h"
#include <functional>
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
  FLUDS(unsigned int num_groups, size_t num_angles, const SPDS& spds)
    : num_groups_(num_groups),
      num_angles_(num_angles),
      num_groups_and_angles_(num_groups_ * num_angles_),
      spds_(spds) {};

  const SPDS& GetSPDS() const { return spds_; }

  virtual void ClearLocalAndReceivePsi() {}
  virtual void ClearSendPsi() {}
  virtual void AllocateInternalLocalPsi() {}
  virtual void AllocateOutgoingPsi() {}

  virtual void AllocateDelayedLocalPsi() {}
  virtual void AllocatePrelocIOutgoingPsi() {}
  virtual void AllocateDelayedPrelocIOutgoingPsi() {}

  std::span<double>& DelayedLocalPsi() { return delayed_local_psi_view_; }
  std::span<double>& DelayedLocalPsiOld() { return delayed_local_psi_old_view_; }
  std::span<double>& DelayedLocalPsiE() { return delayed_local_psiE_view_; }
  std::span<double>& DelayedLocalPsiEOld() { return delayed_local_psiE_old_view_; }
  virtual void SetDelayedLocalPsiOldToNew() {}
  virtual void SetDelayedLocalPsiNewToOld() {}

  std::vector<std::span<double>>& DeplocIOutgoingPsi() { return deplocI_outgoing_psi_view_; }
  std::vector<std::span<double>>& DeplocIOutgoingPsiE() { return deplocI_outgoing_psiE_view_; }

  std::vector<std::span<double>>& PrelocIOutgoingPsi() { return prelocI_outgoing_psi_view_; }
  std::vector<std::span<double>>& PrelocIOutgoingPsiE() { return prelocI_outgoing_psiE_view_; }

  std::vector<std::span<double>>& DelayedPrelocIOutgoingPsi()
  {
    return delayed_prelocI_outgoing_psi_view_;
  }
  std::vector<std::span<double>>& DelayedPrelocIOutgoingPsiE()
  {
    return delayed_prelocI_outgoing_psiE_view_;
  }
  std::vector<std::span<double>>& DelayedPrelocIOutgoingPsiOld()
  {
    return delayed_prelocI_outgoing_psi_old_view_;
  }
  std::vector<std::span<double>>& DelayedPrelocIOutgoingPsiEOld()
  {
    return delayed_prelocI_outgoing_psiE_old_view_;
  }
  virtual void SetDelayedOutgoingPsiOldToNew() {}
  virtual void SetDelayedOutgoingPsiNewToOld() {}

  /// Isotropic nodal correction: (cell local id, face, face node, group index within the
  /// groupset) -> value added to the angular flux in every direction.
  using DelayedFluxCorrection =
    std::function<double(std::uint32_t, unsigned int, unsigned int, unsigned int)>;

  /// Whether this FLUDS stores delayed (lagged) angular fluxes, within or between locations.
  bool HasDelayedAngularFlux() const
  {
    if (not delayed_local_psi_view_.empty())
      return true;
    for (const auto& view : delayed_prelocI_outgoing_psi_view_)
      if (not view.empty())
        return true;
    return false;
  }

  /**
   * Adds an isotropic nodal correction to every new (latest-iterate) delayed angular flux, that is,
   * to the lagged values the next sweep reads where a cyclic dependency was broken, within this
   * location or between locations. The correction is evaluated at the face node of the downwind
   * cell that reads the value.
   *
   * A FLUDS type that stores delayed angular fluxes must override this. The default accepts only a
   * FLUDS without them and otherwise throws, so that lagged fluxes of a type that does not
   * implement the correction fail loudly instead of keeping an error that re-enters on the next
   * sweep.
   */
  virtual void CorrectNewDelayedAngularFlux(const DelayedFluxCorrection& /*correction*/)
  {
    OpenSnLogicalErrorIf(HasDelayedAngularFlux(),
                         "This FLUDS type stores delayed (lagged) angular fluxes but does not "
                         "implement CorrectNewDelayedAngularFlux, so synthetic acceleration cannot "
                         "correct them.");
  }

  virtual ~FLUDS() = default;

protected:
  const unsigned int num_groups_;
  const size_t num_angles_;
  const size_t num_groups_and_angles_;
  const SPDS& spds_;

  std::span<double> delayed_local_psi_view_;
  std::span<double> delayed_local_psi_old_view_;
  std::span<double> delayed_local_psiE_view_;
  std::span<double> delayed_local_psiE_old_view_;
  std::vector<std::span<double>> deplocI_outgoing_psi_view_;
  std::vector<std::span<double>> deplocI_outgoing_psiE_view_;
  std::vector<std::span<double>> prelocI_outgoing_psi_view_;
  std::vector<std::span<double>> prelocI_outgoing_psiE_view_;
  std::vector<std::span<double>> delayed_prelocI_outgoing_psi_view_;
  std::vector<std::span<double>> delayed_prelocI_outgoing_psi_old_view_;
  std::vector<std::span<double>> delayed_prelocI_outgoing_psiE_view_;
  std::vector<std::span<double>> delayed_prelocI_outgoing_psiE_old_view_;
};

} // namespace opensn
