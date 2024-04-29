// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{
namespace lbs
{

class LBSSolver;

void
PowerIterationKEigen(LBSSolver& lbs_solver, double tolerance, int max_iterations, double& k_eff);

} // namespace lbs
} // namespace opensn
