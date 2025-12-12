// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/logging/log_stream.h"
#include "framework/logging/log_exceptions.h"
#include <utility>
#include <vector>
#include <memory>

namespace opensn
{

/**
 * Object for controlling logging.
 *
 * ## Part A: Output logs
 * There are three levels of verbosity in OpenSn: Zero(Default), One and Two.
 * These can be set on the command line via the switch -v followed by a
 * space and the number for the verbosity (0,1 or 2).
 *
 * \code
 * ./opensn InputFile.py -v 1
 * \endcode
 *
 * Printing a log under the auspices of a verbosity level again has
 * numerous options. Firstly, any log can be a normal log, a warning
 * or an error. Secondly, a log can be either location 0 or all the locations
 * in a parallel process environment. The log option enums defined under
 * LOG_LVL are
 *  - LOG_0,                      Used only for location 0
 *  - LOG_0WARNING,               Warning only for location 0
 *  - LOG_0ERROR,                 Error only for location 0
 *  - LOG_0VERBOSE_0,             Default verbosity level
 *  - LOG_0VERBOSE_1,             Used only if verbosity level equals 1
 *  - LOG_0VERBOSE_2,             Used only if verbosity level equals 2
 *  - LOG_ALL,                    Verbose level 0 all locations
 *  - LOG_ALLWARNING,             Warning for any location
 *  - LOG_ALLERROR,               Error for any location
 *  - LOG_ALLVERBOSE_0,     Default verbosity level
 *  - LOG_ALLVERBOSE_1,     Used only if verbosity level equals 1
 *  - LOG_ALLVERBOSE_2,     Used only if verbosity level equals 2
 *
 * A log can be made by first connecting code with the logger. This is done
 * by including the log header and then defining an extern reference to the
 * global object.
 *
 * \code
 * #include "framework/logging/log.h"
 * extern Logger& opensn::log;
 * \endcode
 *
 * A log is then written inside a piece of code as follows:
 *
 * \code
 * void PrintSomethingToLog()
 * {
 *   opensn::log.Log() << "This is printed on location 0 only";
 *   opensn::log.Log0Warning() << "This is a warning";
 *   opensn::log.Log0Error() << "This is an error";
 * }
 * \endcode
 *
 * \verbatim
 * [0]  This is printed on location 0 only
 * [0]  **** WARNING **** This is a warning
 * [0]  ****  ERROR  **** This is an error
 * \endverbatim
 */
class Logger
{
public:
  /// Logging level
  enum LOG_LVL
  {
    LOG_0 = 1,             ///< Used only for location 0
    LOG_0WARNING = 2,      ///< Warning only for location 0
    LOG_0ERROR = 3,        ///< Error only for location 0
    LOG_0VERBOSE_0 = 4,    ///< Default verbosity level
    LOG_0VERBOSE_1 = 5,    ///< Used only if verbosity level equals 1
    LOG_0VERBOSE_2 = 6,    ///< Used only if verbosity level equals 2
    LOG_ALL = 7,           ///< Verbose level 0 all locations
    LOG_ALLWARNING = 8,    ///< Warning for any location
    LOG_ALLERROR = 9,      ///< Error for any location
    LOG_ALLVERBOSE_0 = 10, ///< Default verbosity level
    LOG_ALLVERBOSE_1 = 11, ///< Used only if verbosity level equals 1
    LOG_ALLVERBOSE_2 = 12  ///< Used only if verbosity level equals 2
  };

  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;

  void SetVerbosity(unsigned int level) { verbosity_ = std::min(level, 2U); }
  unsigned int GetVerbosity() const { return verbosity_; }
  LogStream Log(LOG_LVL level = LOG_0);
  LogStream Log0() { return Log(LOG_0); }
  LogStream Log0Warning() { return Log(LOG_0WARNING); }
  LogStream Log0Error() { return Log(LOG_0ERROR); }
  LogStream Log0Verbose0() { return Log(LOG_0VERBOSE_0); }
  LogStream Log0Verbose1() { return Log(LOG_0VERBOSE_1); }
  LogStream Log0Verbose2() { return Log(LOG_0VERBOSE_2); }
  LogStream LogAll() { return Log(LOG_ALL); }
  LogStream LogAllWarning() { return Log(LOG_ALLWARNING); }
  LogStream LogAllError() { return Log(LOG_ALLERROR); }
  LogStream LogAllVerbose0() { return Log(LOG_ALLVERBOSE_0); }
  LogStream LogAllVerbose1() { return Log(LOG_ALLVERBOSE_1); }
  LogStream LogAllVerbose2() { return Log(LOG_ALLVERBOSE_2); }

private:
  Logger() = default;

  DummyStream dummy_stream_;
  unsigned int verbosity_{0};

public:
  static Logger& GetInstance() noexcept
  {
    static Logger instance;
    return instance;
  }
};

} // namespace opensn
