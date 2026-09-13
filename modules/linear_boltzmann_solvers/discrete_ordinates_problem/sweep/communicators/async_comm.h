// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/logging/log.h"
#include <vector>
#include <cstddef>
#include <cstdint>

namespace opensn
{

class SweepCommunicator;
class FLUDS;

class AsynchronousCommunicator
{
public:
  explicit AsynchronousCommunicator(FLUDS& fluds, const SweepCommunicator& sweep_communicator)
    : fluds_(fluds), sweep_communicator_(sweep_communicator)
  {
  }

  virtual ~AsynchronousCommunicator() = default;

protected:
  FLUDS& fluds_;
  const SweepCommunicator& sweep_communicator_;
};

} // namespace opensn
