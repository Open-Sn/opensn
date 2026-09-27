// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>

namespace opensn
{

enum class SteppingMethod
{
  NONE = 0,
  EXPLICIT_EULER = 1,
  IMPLICIT_EULER = 2,
  CRANK_NICOLSON = 3,
  THETA_SCHEME = 4,
};

/// Returns the string name of a time stepping method.
std::string SteppingMethodStringName(SteppingMethod method);
SteppingMethod SteppingMethodFromString(const std::string& name);

/**
 * Returns true if `time` lies in the closed window [start_time, end_time].
 *
 * The comparison allows a relative round-off tolerance so that times accumulated
 * by repeated time steps (e.g. 20 steps of 0.05 giving 1.0000000000000002) are
 * classified the same as the exact time. The tolerance is 1e-12 times the larger magnitude of
 * the evaluation time and the finite bound, with no minimum absolute tolerance. Infinite window
 * bounds are supported.
 */
bool IsTimeInWindow(double time, double start_time, double end_time);

} // namespace opensn
