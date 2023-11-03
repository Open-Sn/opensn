#pragma once

#include <petscksp.h>

namespace lbs
{

PetscErrorCode
GSConvergenceTest(KSP ksp, PetscInt n, PetscReal rnorm, KSPConvergedReason* convergedReason, void*);

} // namespace lbs


