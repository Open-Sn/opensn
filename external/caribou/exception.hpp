/*
 * Created on Tue, April 01 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <sstream>    // std::ostringstream
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string

#include "api_mapping.hpp"  // GPU_API

namespace caribou {

namespace impl {
class CaribouException;
}  // namespace impl

/** @brief Exception to be thrwon when encountering a CUDA/HIP error.*/
class impl::CaribouException : public std::runtime_error {
  public:
    /** @brief Constructor.*/
    explicit CaribouException(::GPU_API(Error_t) error) : std::runtime_error(impl::CaribouException::build_message(error)) {}

  private:
    /** @brief Helper function to build an informative error message.*/
    static std::string build_message(::GPU_API(Error_t) error) {
        std::ostringstream message;
        message << "CUDA/HIP Error(" << static_cast<int>(error) << "): " << ::GPU_API(GetErrorString)(error);
        return message.str();
    }
};

/** @brief Check error and throw.*/
inline void check_error(::GPU_API(Error_t) error) {
    if (error != ::GPU_API(Success)) {
        throw impl::CaribouException(error);
    }
}

}  // namespace caribou
