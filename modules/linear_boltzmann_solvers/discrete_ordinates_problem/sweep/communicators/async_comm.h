// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/logging/log.h"
#include <vector>
#include <cstddef>
#include <cstdint>

namespace opensn
{

class MPICommunicatorSet;
class FLUDS;

class AsynchronousCommunicator
{
public:
  explicit AsynchronousCommunicator(FLUDS& fluds, const MPICommunicatorSet& comm_set)
    : fluds_(fluds), comm_set_(comm_set)
  {
  }

  virtual ~AsynchronousCommunicator() = default;

protected:
  FLUDS& fluds_;
  const MPICommunicatorSet& comm_set_;
};

} // namespace opensn
