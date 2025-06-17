// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <petscksp.h>

namespace opensn
{

PetscErrorCode
GSConvergenceTest(KSP ksp, PetscInt n, PetscReal rnorm, KSPConvergedReason* convergedReason, void*);

} // namespace opensn
