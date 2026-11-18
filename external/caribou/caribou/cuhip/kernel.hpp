/*
 * Created on Mon, November 16 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t
#include <utility>  // std::move

#include "caribou/backend.hpp"
#include "caribou/dim3.hpp"    // caribou::Dim3
#include "caribou/index.hpp"   // caribou::Index
#include "caribou/stream.hpp"  // caribou::Stream

namespace caribou::impl {

// Kernel launcher implementation
// ------------------------------

template <auto Kernel>
class KernelLauncher {
  public:
    KernelLauncher(std::size_t dynamic_shared_mem, caribou::Stream stream) :
    dynamic_shared_mem_(dynamic_shared_mem), stream_(std::move(stream)) {}

    template <typename... Args>
    void launch(const Args &... args) const {
        Kernel<<<this->grid_, this->block_, this->dynamic_shared_mem_, this->stream_>>>(caribou::Index(), args...);
    }

    void set_launch_range(const Dim3 & grid, const Dim3 & block) {
        this->grid_ = static_cast<::dim3>(grid);
        this->block_ = static_cast<::dim3>(block);
    }

  protected:
    ::dim3 grid_;
    ::dim3 block_;
    std::size_t dynamic_shared_mem_;
    mutable caribou::Stream stream_;
};

}  // namespace caribou::impl
