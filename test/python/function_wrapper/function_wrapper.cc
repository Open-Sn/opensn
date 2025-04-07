// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "test/python/src/bindings.h"

using namespace opensn;

namespace unit_tests
{

void
TestCFunction()
{
  opensn::log.Log() << "Hello from a C function";
}

} //  namespace unit_tests
