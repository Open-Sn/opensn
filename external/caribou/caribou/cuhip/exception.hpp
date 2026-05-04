/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <sstream>    // std::ostringstream
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"  // GPU_API

namespace caribou::impl {

class CudaHipException : public std::runtime_error {
  public:
    explicit CudaHipException(::GPU_API(Error_t) error) : std::runtime_error(CudaHipException::build_message(error)) {}

  private:
    static inline std::string build_message(::GPU_API(Error_t) error) {
        std::ostringstream message;
        message << "CUDA/HIP Error(" << static_cast<int>(error) << "): " << ::GPU_API(GetErrorString)(error);
        return message.str();
    }
};

inline void check_error(::GPU_API(Error_t) error) {
    if (error != ::GPU_API(Success)) {
        throw CudaHipException(error);
    }
}

}  // namespace caribou::impl
