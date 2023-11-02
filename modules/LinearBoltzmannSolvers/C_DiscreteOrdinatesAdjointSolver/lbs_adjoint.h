#ifndef LBSADJOINTSOLVER_LBSADJOINT_H
#define LBSADJOINTSOLVER_LBSADJOINT_H

#include <array>

namespace lbs
{
void TestFunction();

std::array<double, 2> MakeExpRepFromP1(const std::array<double, 4>& P1_moments,
                                       bool verbose = false);
} // namespace lbs

#endif // LBSADJOINTSOLVER_LBSADJOINT_H
