#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"

namespace chi_math
{
int LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action);
} // namespace chi_math
