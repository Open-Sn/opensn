#include "linear_matrix_action_Ax.h"

#include "linear_solver_context.h"

#include <petscksp.h>

namespace chi_math
{

template<>
int LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action)
{
  LinearSolverContext<Mat,Vec>* context;
  MatShellGetContext(matrix,&context);

  context->MatrixAction(matrix, vector, //inputs
                        action);        //output

  return 0;
}

}//namespace chi_math

