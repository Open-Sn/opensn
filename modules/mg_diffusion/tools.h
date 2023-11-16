#pragma once

#include "petscksp.h"

namespace opensn
{
namespace mg_diffusion
{
/**My personal monitor monitor.*/
PetscErrorCode MGKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*);
} // namespace mg_diffusion
} // namespace opensn
