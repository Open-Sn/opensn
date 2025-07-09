// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

class LBSProblem;

/// Compute the steady state delayed neutron precursor concentrations.
void ComputePrecursors(LBSProblem& lbs_problem);

} // namespace opensn
