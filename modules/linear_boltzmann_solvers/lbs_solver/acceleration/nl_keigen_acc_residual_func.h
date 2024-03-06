#include <petscsnes.h>

namespace opensn
{
namespace lbs
{

PetscErrorCode NLKEigenAccResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace lbs
} // namespace opensn
