/*
 * Created on Mon, November 16 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "caribou/backend.hpp"
#include BACKEND(index.hpp)  // caribou::impl::Index

namespace caribou {

/**
 * @brief Portable GPU execution index.
 * @details ``Index`` exposes a CUDA-style indexing API for CUDA, HIP and SYCL kernels. Across all backends, axis ``0``
 * is the fastest-varying dimension (CUDA ``x``), axis ``1`` is the middle dimension (CUDA ``y``), and axis ``2`` is the
 * slowest-varying dimension (CUDA ``z``).
 *
 * caribou::KernelLauncher passes an ``Index`` as the first argument of every kernel it launches.
 */
struct Index {
    /// @name Local (thread index in block)
    /// @{
    /**
     * @brief Get the thread index within the current block/work-group along axis ``a``.
     * @tparam a Axis selector: ``0`` fastest, ``1`` middle, ``2`` slowest.
     */
    template <std::uint32_t a>
    requires(a < 3)
    __CRB_DEVICE_FUNC__ std::uint32_t get_local_id(void) const {
        return impl::get_local_id<a>(this->value);
    }
    /**
     * @brief Get the flattened thread index within the current block/work-group.
     * @details The flattening order is CUDA-style, with axis ``0`` fastest and axis ``2`` slowest.
     */
    __CRB_DEVICE_FUNC__ std::uint32_t get_local_linear_id(void) const { return impl::get_local_linear_id(this->value); }
    /**
     * @brief Get the block/work-group extent along axis ``a``.
     * @tparam a Axis selector: ``0`` fastest, ``1`` middle, ``2`` slowest.
     */
    template <std::uint32_t a>
    requires(a < 3)
    __CRB_DEVICE_FUNC__ std::uint32_t get_local_range(void) const {
        return impl::get_local_range<a>(this->value);
    }
    /// @}

    /// @name Group (block index in grid)
    /// @{
    /**
     * @brief Get the block/work-group index within the launch grid along axis ``a``.
     * @tparam a Axis selector: ``0`` fastest, ``1`` middle, ``2`` slowest.
     */
    template <std::uint32_t a>
    requires(a < 3)
    __CRB_DEVICE_FUNC__ std::uint32_t get_group(void) const {
        return impl::get_group<a>(this->value);
    }
    /**
     * @brief Get the flattened block/work-group index within the launch grid.
     * @details The flattening order is CUDA-style, with axis ``0`` fastest and axis ``2`` slowest.
     */
    __CRB_DEVICE_FUNC__ std::uint32_t get_group_linear_id(void) const { return impl::get_group_linear_id(this->value); }
    /**
     * @brief Get the launch-grid extent in blocks/work-groups along axis ``a``.
     * @tparam a Axis selector: ``0`` fastest, ``1`` middle, ``2`` slowest.
     */
    template <std::uint32_t a>
    requires(a < 3)
    __CRB_DEVICE_FUNC__ std::uint32_t get_group_range(void) const {
        return impl::get_group_range<a>(this->value);
    }
    /// @}

    /// @name Global (thread index in grid)
    /// @{
    /**
     * @brief Get the global thread/work-item index along axis ``a``.
     * @tparam a Axis selector: ``0`` fastest, ``1`` middle, ``2`` slowest.
     */
    template <std::uint32_t a>
    requires(a < 3)
    __CRB_DEVICE_FUNC__ std::uint32_t get_global_id(void) const {
        return impl::get_global_id<a>(this->value);
    }
    /**
     * @brief Get the flattened global thread/work-item index.
     * @details The flattening order is CUDA-style, with axis ``0`` fastest and axis ``2`` slowest.
     */
    __CRB_DEVICE_FUNC__ std::uint32_t get_global_linear_id(void) const {
        return impl::get_global_linear_id(this->value);
    }
    /**
     * @brief Get the total launched extent along axis ``a`` (number of groups times the group extent).
     * @tparam a Axis selector: ``0`` fastest, ``1`` middle, ``2`` slowest.
     */
    template <std::uint32_t a>
    requires(a < 3)
    __CRB_DEVICE_FUNC__ std::uint32_t get_global_range(void) const {
        return impl::get_global_range<a>(this->value);
    }
    /// @}

    /// @name Barrier and shared memory
    /// @{
    /**
     * @brief Synchronize all threads of the block/work-group.
     * @tparam WithGlobalMemFence Also make global memory writes visible device-wide before the barrier.
     */
    template <bool WithGlobalMemFence = false>
    __CRB_DEVICE_FUNC__ void sync_block(void) const {
        impl::sync_block<WithGlobalMemFence>(this->value);
    }
    /** @brief Get pointer to the dynamic shared (local) memory configured by caribou::KernelLauncher.*/
    __CRB_DEVICE_FUNC__ char * get_shared_mem(void) const { return this->value.get_shared_mem(); }
    /// @}

    /** @brief Backend-native index payload.*/
    impl::Index value;
};

}  // namespace caribou
