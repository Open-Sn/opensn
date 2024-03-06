#pragma once

#include <petscsnes.h>

namespace opensn
{
namespace lbs
{
PetscErrorCode KEigenSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*);
PetscErrorCode KEigenKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*);
} // namespace lbs
} // namespace opensn
