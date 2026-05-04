// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/logging/log.h"
#include "framework/logging/stringstream_color.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <sstream>

namespace opensn
{

Logger& log = Logger::GetInstance();

LogStream
Logger::Log(LOG_LVL level)
{
  int rank = opensn::mpi_comm.rank();
  std::string header;
  std::ostream* stream = nullptr;
  bool use_dummy = false;
  bool use_color = false;
  auto color = [this](StringStreamColorCode code)
  { return color_enabled_ ? StringStreamColor(code) : std::string{}; };

  switch (level)
  {
    case LOG_0:
    case LOG_0VERBOSE_0:
    {
      if (rank == 0)
      {
        header = "[" + std::to_string(rank) + "]  ";
        stream = &std::cout;
      }
      else
        use_dummy = true;
      break;
    }
    case LOG_0WARNING:
    {
      if (rank == 0)
      {
        header = "[" + std::to_string(rank) + "]  " + color(FG_YELLOW) + "*** WARNING ***  ";
        use_color = color_enabled_;
        stream = &std::cout;
      }
      else
        use_dummy = true;
      break;
    }
    case LOG_0ERROR:
    {
      if (rank == 0)
      {
        header = "[" + std::to_string(rank) + "]  " + color(FG_RED) + "**** ERROR ****  ";
        use_color = color_enabled_;
        stream = &std::cerr;
      }
      else
        use_dummy = true;
      break;
    }
    case LOG_0VERBOSE_1:
    {
      if (rank == 0 and verbosity_ >= 1)
      {
        header = "[" + std::to_string(rank) + "]  " + color(FG_CYAN);
        use_color = color_enabled_;
        stream = &std::cout;
      }
      else
        use_dummy = true;
      break;
    }
    case LOG_0VERBOSE_2:
    {
      if (rank == 0 and verbosity_ >= 2)
      {
        header = "[" + std::to_string(rank) + "]  " + color(FG_MAGENTA);
        use_color = color_enabled_;
        stream = &std::cout;
      }
      else
        use_dummy = true;
      break;
    }
    case LOG_ALL:
    case LOG_ALLVERBOSE_0:
    {
      header = "[" + std::to_string(rank) + "]  ";
      stream = &std::cout;
      break;
    }
    case LOG_ALLWARNING:
    {
      header = "[" + std::to_string(rank) + "]  " + color(FG_YELLOW) + "*** WARNING ***  ";
      use_color = color_enabled_;
      stream = &std::cout;
      break;
    }
    case LOG_ALLERROR:
    {
      header = "[" + std::to_string(rank) + "]  " + color(FG_RED) + "**** ERROR ****  ";
      use_color = color_enabled_;
      stream = &std::cerr;
      break;
    }
    case LOG_ALLVERBOSE_1:
    {
      if (verbosity_ >= 1)
      {
        header = "[" + std::to_string(rank) + "]  " + color(FG_CYAN);
        use_color = color_enabled_;
        stream = &std::cout;
      }
      else
        use_dummy = true;
      break;
    }
    case LOG_ALLVERBOSE_2:
    {
      if (verbosity_ >= 2)
      {
        header = "[" + std::to_string(rank) + "]  " + color(FG_MAGENTA);
        use_color = color_enabled_;
        stream = &std::cout;
      }
      else
        use_dummy = true;
      break;
    }
    default:
      use_dummy = true;
      break;
  }

  if (use_dummy)
    return {&dummy_stream_, " ", true};

  return {stream, header, false, use_color};
}

} // namespace opensn
