#pragma once

#include <petscsnes.h>

namespace lbs
{

PetscErrorCode NLKEigenResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace lbs
