#ifndef CHITECH_TOOLS_H
#define CHITECH_TOOLS_H

#include "petscksp.h"

namespace mg_diffusion
{
PetscErrorCode MGKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*);
}

#endif // CHITECH_TOOLS_H
