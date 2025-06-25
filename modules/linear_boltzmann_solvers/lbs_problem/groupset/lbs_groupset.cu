// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"

namespace opensn
{

void
LBSGroupset::InitializeGPUCarriers()
{
  QuadratureCarrier* quad = new QuadratureCarrier(*this);
  quad_carrier = quad;
}

void
LBSGroupset::ResetGPUCarriers()
{
  if (quad_carrier)
  {
    QuadratureCarrier* quad = reinterpret_cast<QuadratureCarrier*>(quad_carrier);
    delete quad;
    quad_carrier = nullptr;
  }
}

} // namespace opensn
