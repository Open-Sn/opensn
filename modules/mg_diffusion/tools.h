// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "petscksp.h"

namespace opensn
{
namespace mg_diffusion
{
/**My personal monitor monitor.*/
PetscErrorCode MGKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*);
} // namespace mg_diffusion
} // namespace opensn
