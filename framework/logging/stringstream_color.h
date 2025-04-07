// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <unistd.h>

namespace opensn
{

enum StringStreamColorCode
{
  RESET = 0,
  FG_BOLD = 1,
  FG_UNDERLINE = 4,
  FG_BOLD_OFF = 21,
  FG_UNDERLINE_OFF = 24,
  FG_RED = 31,
  FG_GREEN = 32,
  FG_YELLOW = 33,
  FG_BLUE = 34,
  FG_MAGENTA = 35,
  FG_CYAN = 36,
  FG_WHITE = 37,
  FG_DEFAULT = 39,
};

std::string StringStreamColor(StringStreamColorCode code);

} // namespace opensn
