// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <petscsnes.h>

namespace opensn
{

/**This function evaluates the flux moments based k-eigenvalue transport
 * residual of the form
 * \f$ r(\phi) = DL^{-1} (\frac{1}{k} F\phi + MS \phi) - \phi \f$.
 */
PetscErrorCode NLKEigenResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace opensn
