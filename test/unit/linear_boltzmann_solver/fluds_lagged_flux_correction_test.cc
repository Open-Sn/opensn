// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{

/// A FLUDS type with optional lagged storage that does not implement the lagged-flux correction.
class UncorrectedFLUDS : public FLUDS
{
public:
  UncorrectedFLUDS(const SPDS& spds, bool lagged_local, bool lagged_nonlocal)
    : FLUDS(1, 1, spds), local_(2, 0.0), nonlocal_(3, 0.0)
  {
    if (lagged_local)
      delayed_local_psi_view_ = local_;
    // One upwind location, with or without lagged values.
    if (lagged_nonlocal)
      delayed_prelocI_outgoing_psi_view_.emplace_back(nonlocal_);
    else
      delayed_prelocI_outgoing_psi_view_.emplace_back();
  }

private:
  std::vector<double> local_;
  std::vector<double> nonlocal_;
};

double
UnitCorrection(std::uint32_t /*cell_local_id*/,
               unsigned int /*face*/,
               unsigned int /*face_node*/,
               unsigned int /*group_idx*/)
{
  return 1.0;
}

TEST(FLUDSLaggedFluxCorrection, WithoutLaggedFluxIsNoOp)
{
  const SPDS spds(Vector3(0.0, 0.0, 1.0), nullptr);
  UncorrectedFLUDS fluds(spds, false, false);
  EXPECT_FALSE(fluds.HasDelayedAngularFlux());
  EXPECT_NO_THROW(fluds.CorrectNewDelayedAngularFlux(UnitCorrection));
}

TEST(FLUDSLaggedFluxCorrection, LocalLaggedFluxWithoutCorrectionThrows)
{
  const SPDS spds(Vector3(0.0, 0.0, 1.0), nullptr);
  UncorrectedFLUDS fluds(spds, true, false);
  EXPECT_TRUE(fluds.HasDelayedAngularFlux());
  EXPECT_THROW(fluds.CorrectNewDelayedAngularFlux(UnitCorrection), std::logic_error);
}

TEST(FLUDSLaggedFluxCorrection, NonlocalLaggedFluxWithoutCorrectionThrows)
{
  const SPDS spds(Vector3(0.0, 0.0, 1.0), nullptr);
  UncorrectedFLUDS fluds(spds, false, true);
  EXPECT_TRUE(fluds.HasDelayedAngularFlux());
  EXPECT_THROW(fluds.CorrectNewDelayedAngularFlux(UnitCorrection), std::logic_error);
}

} // namespace
