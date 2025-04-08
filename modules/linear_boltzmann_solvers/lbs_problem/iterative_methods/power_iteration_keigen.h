// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

class LBSProblem;

void PowerIterationKEigenSolver(LBSProblem& lbs_problem,
                                double tolerance,
                                int max_iterations,
                                double& k_eff);

} // namespace opensn
