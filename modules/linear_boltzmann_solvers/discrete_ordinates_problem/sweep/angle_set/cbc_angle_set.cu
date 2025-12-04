// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/storage.h"
#include "caribou/cuda/stream.hpp"
#include "caribou/caribou.h"

namespace opensn
{

void
CBC_AngleSet::AssociateAngleSetWithDeviceStructures()
{
  stream_ = caribou::Stream{};
  CBCD_FLUDS* cbcd_fluds = dynamic_cast<CBCD_FLUDS*>(fluds_.get());
  cbcd_fluds->SetAngleSet(*this);
}

} // namespace opensn