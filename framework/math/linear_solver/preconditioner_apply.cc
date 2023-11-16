#include "framework/math/linear_solver/preconditioner_apply.h"

#include "framework/math/linear_solver/preconditioner_context.h"

#include <petscksp.h>

namespace opensn
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

} // namespace opensn
