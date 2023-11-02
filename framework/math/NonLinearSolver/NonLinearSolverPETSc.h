#ifndef CHITECH_NONLINEARSOLVERPETSC_H
#define CHITECH_NONLINEARSOLVERPETSC_H

#include "NonLinearSolver.h"

#include <petscsnes.h>

namespace chi_math
{

using NonLinearSolverPETSc = NonLinearSolver<Mat, Vec, SNES>;

}

#endif // CHITECH_NONLINEARSOLVERPETSC_H
