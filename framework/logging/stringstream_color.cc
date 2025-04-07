// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/logging/stringstream_color.h"
#include "framework/runtime.h"

namespace opensn
{

std::string
StringStreamColor(StringStreamColorCode code)
{
  if (suppress_color)
    return {};
  return "\033[" + std::to_string(static_cast<int>(code)) + "m";
}

} // namespace opensn
