// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include <petscsnes.h>

namespace opensn
{

PetscErrorCode NLKEigenAccResidualFunction(SNES snes, Vec phi, Vec r, void* ctx);

} // namespace opensn
