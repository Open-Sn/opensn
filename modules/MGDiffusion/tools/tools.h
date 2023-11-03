#pragma once

#include "petscksp.h"

namespace mg_diffusion
{
PetscErrorCode MGKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*);
}
