// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>

namespace opensn
{

struct KResidualFunctionContext
{
  std::string solver_name;
  double k_eff = 1.0;
};

} // namespace opensn
