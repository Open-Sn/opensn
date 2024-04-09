// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/logging/log_stream.h"
#include "framework/logging/stringstream_color.h"

namespace opensn
{

LogStream::~LogStream()
{
  if (dummy_)
    return;

  std::string line, oline;
  while (std::getline(*this, line))
    oline += log_header_ + line + '\n' + StringStreamColor(RESET);

  if (not oline.empty())
    *log_stream_ << oline << std::flush;
}

} // namespace opensn
