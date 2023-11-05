#include "opensn/framework/math/LinearSolver/preconditioner_apply.h"

#include "opensn/framework/math/LinearSolver/preconditioner_context.h"

#include <petscksp.h>

namespace chi_math
{

template <>
int
PreconditionerApplication(PC pc, Vec vector, Vec action)
{
  PreconditionerContext<PC, Vec>* context;
  PCShellGetContext(pc, &context);

  context->PCApply(pc, vector, action);

  return 0;
}

} // namespace chi_math
