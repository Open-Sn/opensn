// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"

namespace opensn
{

void
LBSGroupset::InitializeQuadratureCarrier()
{
  quad_carrier = std::make_shared<QuadratureCarrier>(*this);
}

void
LBSGroupset::ResetQuadratureCarrier()
{
  quad_carrier.reset();
}

} // namespace opensn
