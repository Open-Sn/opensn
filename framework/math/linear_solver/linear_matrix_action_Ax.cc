#include "framework/math/linear_solver/linear_matrix_action_Ax.h"

#include <petscksp.h>

namespace opensn
{

int
LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action)
{
  LinearSolverContext* context;
  MatShellGetContext(matrix, &context);

  context->MatrixAction(matrix, vector, action);

  return 0;
}

} // namespace opensn
