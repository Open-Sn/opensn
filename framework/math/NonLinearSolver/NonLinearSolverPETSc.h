#pragma once

#include "NonLinearSolver.h"

#include <petscsnes.h>

namespace chi_math
{

using NonLinearSolverPETSc = NonLinearSolver<Mat, Vec, SNES>;

}
