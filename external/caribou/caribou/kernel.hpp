/*
 * Created on Mon, November 16 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstddef>      // std::size_t
#include <type_traits>  // std::is_trivially_copyable_v, std::remove_cvref_t
#include <utility>      // std::move

#include "caribou/backend.hpp"
#include "caribou/dim3.hpp"    // caribou::Dim3
#include "caribou/index.hpp"   // caribou::Index
#include "caribou/stream.hpp"  // caribou::Stream
#include BACKEND(kernel.hpp)   // caribou::impl::KernelLauncher

namespace caribou {

/**
 * @brief Backend-agnostic launcher for a compile-time GPU kernel.
 * @details This wrapper exposes a CUDA/HIP-like launch configuration while dispatching to the active GPU backend. For
 * CUDA/HIP, ``Kernel`` is launched as a global kernel with the configured grid, block, dynamic shared memory and
 * stream. For SYCL, ``Kernel`` is invoked from a ``parallel_for`` over the equivalent ``sycl::nd_range<3>``. In all
 * backends, the launcher injects a caribou::Index value as the first kernel argument, so kernels can query local,
 * group and global indices through a unified API.
 *
 * Example:
 * @code {.cu}
 * __CRB_GLOBAL_FUNC__ void saxpy_kernel(caribou::Index index, std::uint32_t n, double alpha, const double * x,
 *                                       double * y) {
 *     const std::uint32_t i = index.get_global_id<0>();
 *     if (i < n) {
 *         y[i] += alpha * x[i];
 *     }
 * }
 *
 * void saxpy(std::uint32_t n, double alpha, const double * x, double * y, caribou::Stream & stream) {
 *     const caribou::Dim3 block(128);
 *     const caribou::Dim3 grid((n + block.x - 1) / block.x);
 *     caribou::KernelLauncher<saxpy_kernel>(grid, block, 0, stream)(n, alpha, x, y);
 * }
 * @endcode
 * @tparam Kernel Pointer to the kernel entry point. The kernel takes caribou::Index as its first argument, followed by
 * trivially-copyable user arguments.
 */
template <auto Kernel>
class KernelLauncher : public impl::KernelLauncher<Kernel> {
  public:
    /**
     * @brief Construct a launcher with a fixed execution configuration.
     * @param grid Grid dimensions (number of blocks) in CUDA/HIP ordering.
     * @param block Block dimensions in CUDA/HIP ordering.
     * @param dynamic_shared_mem Dynamic shared-memory size in bytes.
     * @param stream Stream or queue used for kernel submission.
     * @note For SYCL, the grid and block dimensions are converted internally to the global and local ranges of a
     * ``sycl::nd_range<3>``.
     */
    KernelLauncher(const Dim3 & grid, const Dim3 & block, std::size_t dynamic_shared_mem = 0,
                   Stream stream = Stream::get_null_stream()) :
    impl::KernelLauncher<Kernel>(dynamic_shared_mem, std::move(stream)) {
        this->set_grid_block(grid, block);
    }

    /**
     * @brief Launch the bound kernel with the stored configuration.
     * @details The launch is asynchronous with respect to the host. A caribou::Index is prepended to @p args.
     * @param args Kernel arguments forwarded to ``Kernel`` (copied by value to the device).
     */
    template <typename... Args>
    requires(std::is_trivially_copyable_v<std::remove_cvref_t<Args>> && ...)
    void operator()(const Args &... args) const {
        this->launch(args...);
    }

    /**
     * @brief Update the launch grid and block dimensions.
     * @param grid Grid dimensions (number of blocks) in CUDA/HIP ordering.
     * @param block Block dimensions in CUDA/HIP ordering.
     */
    void set_grid_block(const Dim3 & grid, const Dim3 & block) { this->set_launch_range(grid, block); }

    /**
     * @brief Set the dynamic shared-memory size used for subsequent launches.
     * @note For SYCL, nonzero dynamic shared memory is mapped to a ``sycl::local_accessor<char, 1>``.
     */
    void set_shared_mem(std::size_t dynamic_shared_mem) { this->dynamic_shared_mem_ = dynamic_shared_mem; }
};

}  // namespace caribou
