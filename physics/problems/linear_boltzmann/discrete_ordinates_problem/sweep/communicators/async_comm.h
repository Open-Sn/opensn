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

  virtual std::vector<double>& InitGetDownwindMessageData(int location_id,
                                                          uint64_t cell_global_id,
                                                          unsigned int face_id,
                                                          size_t angle_set_id,
                                                          size_t data_size)
  {
    OpenSnLogicalError("Method not implemented");
  }

protected:
  FLUDS& fluds_;
  const MPICommunicatorSet& comm_set_;
};

} // namespace opensn
