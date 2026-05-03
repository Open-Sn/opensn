// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/math/functions/function.h"
#include "framework/data_types/ndarray.h"
#include <map>
#include <memory>
#include <vector>

namespace opensn
{

/// Arbitrary incident angular flux specified by a user-supplied function.
///
/// The AngularFluxFunction is expected to return the incoming angular flux for a given
/// energy group index and quadrature direction index. The returned values are uniform for
/// all faces with this boundary id.
class ArbitraryBoundary : public SweepBoundary
{
public:
  ArbitraryBoundary(BoundaryBank& bank,
                    const std::vector<LBSGroupset>& groupsets,
                    std::shared_ptr<AngularFluxTimeFunction> angular_flux_function);

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id,
                      unsigned int group_idx) override;

  void UpdateBoundaryFlux(const std::vector<LBSGroupset>& groupsets) override;

  void InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets) override;
  std::uint64_t
  GetOffsetToAngleset(const FaceNode& face_node, AngleSet& anglset, bool is_outgoing) override;

private:
  /// Function that dictates the angular flux.
  std::shared_ptr<AngularFluxTimeFunction> angular_flux_function_;

  /// Additional data specific to arbitrary boundary for each groupset.
  struct ExtraData
  {
    /// Map from angle number in a quadrature to the internal angle index per groupset.
    std::vector<std::uint64_t> map_dirnum;
  };
  /// List of data per groupset.
  std::vector<ExtraData> extra_data_;
};

} // namespace opensn
