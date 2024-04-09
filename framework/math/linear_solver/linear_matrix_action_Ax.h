// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"

namespace opensn
{

int LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action);

} // namespace opensn
