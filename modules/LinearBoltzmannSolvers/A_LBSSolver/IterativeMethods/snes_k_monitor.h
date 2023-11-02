#ifndef CHITECH_LBS_SNES_MONITOR_H
#define CHITECH_LBS_SNES_MONITOR_H

#include <petscsnes.h>

namespace lbs
{
  PetscErrorCode
  KEigenSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*);
  PetscErrorCode
  KEigenKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *);
}//namespace lbs

#endif //CHITECH_LBS_SNES_MONITOR_H
