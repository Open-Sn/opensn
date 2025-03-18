// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <petscsnes.h>

namespace opensn
{

PetscErrorCode KEigenSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*);
PetscErrorCode KEigenKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*);

} // namespace opensn
