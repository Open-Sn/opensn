// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

class LBSProblem;

/// Runs power iterations on `lbs_problem`, starting from k_eff = 1, until k_eff changes by less
/// than `tolerance` or `max_iterations` is reached, and returns the result in `k_eff`.
void PowerIterationKEigen(LBSProblem& lbs_problem,
                          double tolerance,
                          unsigned int max_iterations,
                          double& k_eff);

} // namespace opensn
