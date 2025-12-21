// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "external/caribou/caribou.h"
#include "external/caribou/cuda/stream.hpp"

namespace opensn
{
void
CBC_AngleSet::CreateStream()
{
  if (stream_ == nullptr)
  {
    caribou::Stream* stream = new caribou::Stream();
    stream_ = stream;
  }
}

void
CBC_AngleSet::DestroyStream()
{
  if (stream_ != nullptr)
  {
    caribou::Stream* stream = static_cast<caribou::Stream*>(stream_);
    delete stream;
    stream_ = nullptr;
  }
}

void
CBC_AngleSet::AssociateAngleSetWithFLUDS()
{
  CBCD_FLUDS* cbcd_fluds = dynamic_cast<CBCD_FLUDS*>(fluds_.get());
  cbcd_fluds->SetAngleSet(*this);
}

} // namespace opensn