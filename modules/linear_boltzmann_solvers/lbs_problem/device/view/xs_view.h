// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/inline_macro.h"
#include <cstdint>

namespace opensn
{

/// Cross section view from contiguous block of memory
struct XSView
{
  /// Get the total cross section view at a given index in the cross section array.
  __inline_host_dev__ XSView(char* xs_data, std::uint32_t index)
  {
    std::uint32_t* num_groups_data = reinterpret_cast<std::uint32_t*>(xs_data) + 1;
    num_groups = *(num_groups_data++);
    total_xs_data = reinterpret_cast<double*>(num_groups_data);
    total_xs_data += index * num_groups;
  }
  double* total_xs_data;
  std::uint32_t num_groups;
};

} // namespace opensn
