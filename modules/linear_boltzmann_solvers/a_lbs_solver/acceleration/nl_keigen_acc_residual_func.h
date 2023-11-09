#include <petscsnes.h>

namespace lbs::acceleration
{

PetscErrorCode NLKEigenAccResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace lbs::acceleration