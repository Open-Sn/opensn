#pragma once

#include <string>
#include <stdexcept>

#define ChiInvalidArgumentIf(condition, message)                                                   \
  if (condition) throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": " + message)
#define ChiInvalidArgument(message)                                                                \
  throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define ChiLogicalErrorIf(condition, message)                                                      \
  if (condition) throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define ChiLogicalError(message)                                                                   \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": " + message)

#define ChiRecoverableInvalidArgument(condition, message)                                          \
  {                                                                                                \
    if (condition)                                                                                 \
      throw std::RecoverableException(std::string("Recoverable Invalid Argument: "),               \
                                      std::string(__PRETTY_FUNCTION__) + ": " + #message);         \
  }

#define ChiRecoverableLogicalError(condition, message)                                             \
  {                                                                                                \
    if (condition)                                                                                 \
      throw std::RecoverableException(std::string("Recoverable Logic Error: ")                     \
                                        std::string(__PRETTY_FUNCTION__) +                         \
                                      ": " + #message);                                            \
  }


