#pragma once

#include "framework/math/nonlinear_solver/nonlinear_solver.h"

#include <petscsnes.h>

namespace chi_math
{

using NonLinearSolverPETSc = NonLinearSolver<Mat, Vec, SNES>;

}
