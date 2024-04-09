// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>

namespace opensn
{
namespace lbs
{
void TestFunction();

std::array<double, 2> MakeExpRepFromP1(const std::array<double, 4>& P1_moments,
                                       bool verbose = false);
} // namespace lbs
} // namespace opensn
