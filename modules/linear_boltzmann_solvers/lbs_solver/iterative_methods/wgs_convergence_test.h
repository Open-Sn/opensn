#pragma once

#include <petscksp.h>

namespace opensn
{
namespace lbs
{

PetscErrorCode
GSConvergenceTest(KSP ksp, PetscInt n, PetscReal rnorm, KSPConvergedReason* convergedReason, void*);

} // namespace lbs
} // namespace opensn
