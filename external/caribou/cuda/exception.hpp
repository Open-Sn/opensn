/*
 * Created on Tue, April 01 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <sstream>    // std::ostringstream
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string

namespace caribou {

namespace impl {
class CudaException;
}  // namespace impl

/** @brief Exception to be thrwon when encountering a CUDA error.*/
class impl::CudaException : public std::runtime_error {
  public:
    /** @brief Constructor.*/
    explicit CudaException(::cudaError_t error) : std::runtime_error(impl::CudaException::build_message(error)) {}

  private:
    /** @brief Helper function to build an informative error message.*/
    static std::string build_message(::cudaError_t error) {
        std::ostringstream message;
        message << "CUDA Error(" << static_cast<int>(error) << "): " << ::cudaGetErrorString(error);
        return message.str();
    }
};

/** @brief Check error and throw.*/
inline void check_error(::cudaError_t error) {
    if (error != 0) {
        throw impl::CudaException(error);
    }
}

}  // namespace caribou
