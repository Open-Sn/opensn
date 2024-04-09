// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/logging/log.h"
#include "framework/logging/stringstream_color.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <sstream>

namespace opensn
{

Logger&
Logger::GetInstance() noexcept
{
  static Logger instance;
  return instance;
}

Logger::Logger() noexcept
{
  verbosity_ = 0;
}

LogStream
Logger::Log(LOG_LVL level)
{
  switch (level)
  {
    case LOG_0VERBOSE_0:
    case LOG_0:
    {
      if (opensn::mpi_comm.rank() == 0)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_0WARNING:
    {
      if (opensn::mpi_comm.rank() == 0)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_YELLOW) + "**WARNING** ";
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_0ERROR:
    {
      if (opensn::mpi_comm.rank() == 0)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_RED) + "**!**ERROR**!** ";
        return {&std::cerr, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }

    case LOG_0VERBOSE_1:
    {
      if ((opensn::mpi_comm.rank() == 0) and (verbosity_ >= 1))
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_CYAN);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case Logger::LOG_LVL::LOG_0VERBOSE_2:
    {
      if ((opensn::mpi_comm.rank() == 0) and (verbosity_ >= 2))
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_MAGENTA);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_ALLVERBOSE_0:
    case LOG_ALL:
    {
      std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
      return {&std::cout, header};
    }
    case LOG_ALLWARNING:
    {
      std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
      header += StringStreamColor(FG_YELLOW) + "**WARNING** ";
      return {&std::cout, header};
    }
    case LOG_ALLERROR:
    {
      std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
      header += StringStreamColor(FG_RED) + "**!**ERROR**!** ";
      return {&std::cerr, header};
    }

    case LOG_ALLVERBOSE_1:
    {
      if (verbosity_ >= 1)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_CYAN);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_ALLVERBOSE_2:
    {
      if (verbosity_ >= 2)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_MAGENTA);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    default:
      std::string header = " ";
      return {&dummy_stream_, header};
  }
}

void
Logger::SetVerbosity(int int_level)
{
  verbosity_ = std::min(int_level, 2);
}

int
Logger::GetVerbosity() const
{
  return verbosity_;
}

} // namespace opensn
