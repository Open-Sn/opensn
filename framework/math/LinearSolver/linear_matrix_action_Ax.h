#ifndef CHITECH_MATRIX_ACTION_AX_H
#define CHITECH_MATRIX_ACTION_AX_H

#include "linear_solver_context.h"

namespace chi_math
{
template<class MatType, class VecType>
int LinearSolverMatrixAction(MatType matrix, VecType vector, VecType action);
}//namespace chi_math

#endif //CHITECH_MATRIX_ACTION_AX_H
