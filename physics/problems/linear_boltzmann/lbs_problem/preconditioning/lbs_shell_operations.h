// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/solvers/iterative_methods/wgs_context.h"
#include <petscksp.h>

namespace opensn
{

/// Applies WGDSA or TGDSA to the given input vector.
int WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);
/// Applies WGDSA or TGDSA to the given input vector.
int WGDSA_TGDSA_PreConditionerMult2(WGSContext& gs_context_ptr, Vec phi_input, Vec pc_output);
/// Applies TGDSA to the given input vector.
int MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);

} // namespace opensn
