// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver.h"

namespace opensnlua
{

void SolverInitialize(std::shared_ptr<opensn::Solver> solver);
void SolverExecute(std::shared_ptr<opensn::Solver> solver);
void SolverStep(std::shared_ptr<opensn::Solver> solver);
void SolverAdvance(std::shared_ptr<opensn::Solver> solver);

} // namespace opensnlua
