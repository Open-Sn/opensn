// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <stdexcept>

#define OpenSnInvalidArgumentIf(condition, message)                                                \
  if (condition)                                                                                   \
  throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define OpenSnInvalidArgument(message)                                                             \
  throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define OpenSnLogicalErrorIf(condition, message)                                                   \
  if (condition)                                                                                   \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define OpenSnLogicalError(message)                                                                \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define OpenSnRecoverableInvalidArgument(condition, message)                                       \
  {                                                                                                \
    if (condition)                                                                                 \
      throw std::RecoverableException(std::string("Recoverable Invalid Argument: "),               \
                                      std::string(__PRETTY_FUNCTION__) + ": " + #message);         \
  }

#define OpenSnRecoverableLogicalError(condition, message)                                          \
  {                                                                                                \
    if (condition)                                                                                 \
      throw std::RecoverableException(std::string("Recoverable Logic Error: ")                     \
                                        std::string(__PRETTY_FUNCTION__) +                         \
                                      ": " + #message);                                            \
  }
