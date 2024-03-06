#pragma once

#include <petscsnes.h>

namespace opensn
{
namespace lbs
{

/**This function evaluates the flux moments based k-eigenvalue transport
 * residual of the form
 * \f$ r(\phi) = DL^{-1} (\frac{1}{k} F\phi + MS \phi) - \phi \f$.
 */
PetscErrorCode NLKEigenResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace lbs
} // namespace opensn
