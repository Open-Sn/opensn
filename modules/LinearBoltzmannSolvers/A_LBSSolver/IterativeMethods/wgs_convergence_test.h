#ifndef CHITECH_WGS_CONVERGENCE_TEST_H
#define CHITECH_WGS_CONVERGENCE_TEST_H

#include <petscksp.h>

namespace lbs
{

PetscErrorCode GSConvergenceTest(
  KSP ksp, PetscInt n, PetscReal rnorm,
  KSPConvergedReason* convergedReason, void*);

}//namespace lbs

#endif //CHITECH_WGS_CONVERGENCE_TEST_H
