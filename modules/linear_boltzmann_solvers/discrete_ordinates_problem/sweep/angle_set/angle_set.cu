// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"

namespace opensn
{

void
AngleSet::InitializeMemoryPin()
{
  MemoryPinner<std::uint32_t>* directions = new MemoryPinner<std::uint32_t>(angles_);
  directions->CopyToDevice();
  memory_pin_ = directions;
}

void
AngleSet::ResetMemoryPin()
{
  if (memory_pin_)
  {
    auto* directions = reinterpret_cast<MemoryPinner<std::uint32_t>*>(memory_pin_);
    delete directions;
    memory_pin_ = nullptr;
  }
}

} // namespace opensn
