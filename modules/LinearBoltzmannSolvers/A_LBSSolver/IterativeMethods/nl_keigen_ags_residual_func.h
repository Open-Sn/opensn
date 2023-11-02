#ifndef CHITECH_LBS_NL_KEIGEN_AGS_RESIDUAL_FUNC_H
#define CHITECH_LBS_NL_KEIGEN_AGS_RESIDUAL_FUNC_H

#include <petscsnes.h>

namespace lbs
{

PetscErrorCode NLKEigenResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace lbs

#endif // CHITECH_LBS_NL_KEIGEN_AGS_RESIDUAL_FUNC_H
