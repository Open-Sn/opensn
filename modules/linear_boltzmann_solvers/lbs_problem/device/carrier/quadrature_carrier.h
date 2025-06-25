// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <cstdint>

namespace opensn
{

/// Object managing angleset and quadrature on GPU.
class QuadratureCarrier : public Carrier
{
public:
  QuadratureCarrier(const LBSGroupset& groupset);

protected:
  std::uint64_t ComputeSize(const LBSGroupset& groupset);
  void Assemble(const LBSGroupset& groupset);
};

} // namespace opensn
