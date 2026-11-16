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
#include "caribou/index.hpp"   // caribou::Index, caribou::impl::Index
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
        if (this->dynamic_shared_mem_ == 0) {
            this->stream_.parallel_for(::sycl::nd_range<3>(this->global_, this->local_),
                                       [=](::sycl::nd_item<3> work_index) {
                                           Kernel(caribou::Index{impl::Index{work_index, nullptr}}, args...);
                                       });
        } else {
            this->stream_.submit([&](::sycl::handler & handler) {
                ::sycl::local_accessor<char, 1> crb_shared_mem(::sycl::range<1>(this->dynamic_shared_mem_), handler);
                handler.parallel_for(::sycl::nd_range<3>(this->global_, this->local_),
                                     [=](::sycl::nd_item<3> work_index) {
                                         Kernel(caribou::Index{impl::Index{work_index, &crb_shared_mem[0]}}, args...);
                                     });
            });
        }
    }

    void set_launch_range(const Dim3 & grid, const Dim3 & block) {
        this->local_ = static_cast<::sycl::range<3>>(block);
        this->global_ = static_cast<::sycl::range<3>>(grid) * this->local_;
    }

  protected:
    ::sycl::range<3> global_{0, 0, 0};
    ::sycl::range<3> local_{0, 0, 0};
    std::size_t dynamic_shared_mem_;
    mutable caribou::Stream stream_;
};

}  // namespace caribou::impl
