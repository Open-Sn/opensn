#pragma once

#include <array>

namespace lbs
{
void TestFunction();

std::array<double, 2> MakeExpRepFromP1(const std::array<double, 4>& P1_moments,
                                       bool verbose = false);
} // namespace lbs
