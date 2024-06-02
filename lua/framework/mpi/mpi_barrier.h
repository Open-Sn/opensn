// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/lua.h"

namespace opensnlua
{

/** Blocks until all processes in the communicator have reached this routine.
 *
 * \ingroup chiMPI
 * \author Jan
 */
int MPIBarrier(lua_State* L);

} // namespace opensnlua
