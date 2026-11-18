/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t
#include <string>   // std::string

#include "caribou/backend.hpp"
#include BACKEND(device.hpp)

namespace caribou {

/** @brief Get number of GPUs associated to the current host.*/
inline std::uint32_t get_num_gpus(void) { return impl::get_num_gpus(); }

/**
 * @brief Get the index of the GPU device on which subsequent operations of the calling thread are performed.
 * @return Index of the device in ``[0, caribou::get_num_gpus())``.
 */
inline std::uint32_t get_current_device(void) { return impl::get_current_device(); }

/**
 * @brief Set the GPU device on which subsequent operations of the calling thread are performed.
 * @details Memory allocations, synchronous copies, newly created streams and the null stream returned by
 * caribou::Stream::get_null_stream all target the selected device.
 * @note For SYCL implementation, there is no device binding in the backend. The current device is a per-thread index
 * tracked by this library, selecting the default queue and the target device of the operations submitted afterwards.
 * @param device_id Index of the device to switch to, in ``[0, caribou::get_num_gpus())``.
 */
inline void set_current_device(std::uint32_t device_id) { impl::set_current_device(device_id); }

/**
 * @brief Get the PCI bus ID of GPU @p device_id, which identifies the GPU on its node.
 * @details The ID is formatted as ``DDDD:BB:DD.F`` (domain, bus, device, function).
 * @note For SYCL implementation, the query relies on the ``sycl_ext_intel_device_info`` extension. A
 * ``sycl::exception`` with ``sycl::errc::feature_not_supported`` is thrown if the device does not support it.
 * @param device_id Index of the device, in ``[0, caribou::get_num_gpus())``.
 */
inline std::string get_pci_bus_id(std::uint32_t device_id) { return impl::get_pci_bus_id(device_id); }

/**
 * @brief Get the number of compute units of GPU @p device_id.
 * @details Streaming multiprocessors (SMs) for CUDA, compute units (CUs) for HIP, and
 * ``sycl::info::device::max_compute_units`` for SYCL (Xe cores or EUs on Intel GPUs, SMs on NVIDIA GPUs).
 * @param device_id Index of the device, in ``[0, caribou::get_num_gpus())``.
 */
inline std::uint32_t get_num_compute_units(std::uint32_t device_id) { return impl::get_num_compute_units(device_id); }

/**
 * @brief Get warp/subgroup size.
 * @note For SYCL implementation, this is fixed as a compile-time constant. However, developper
 * should profile kernel by picking one of the values in caribou::impl::get_subgroup_size. There is
 * no universal optimal value.
 */
inline std::uint32_t get_warp_size(void) { return impl::get_warp_size(); }

}  // namespace caribou
