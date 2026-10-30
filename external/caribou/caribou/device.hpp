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

/**
 * @brief Get number of GPUs that the current process can use.
 * @note Returns 0 when the runtime finds no GPU or no usable driver.
 */
inline std::uint32_t get_num_gpus(void) { return impl::get_num_gpus(); }

/**
 * @brief Make GPU @p device_id the current device of the calling host thread.
 * @note Not supported by the SYCL backend, which selects its device when the library loads.
 */
inline void set_device(std::uint32_t device_id) { impl::set_device(device_id); }

/**
 * @brief Get the PCI bus ID of GPU @p device_id, which identifies the GPU on its node.
 * @note Not supported by the SYCL backend.
 */
inline std::string get_pci_bus_id(std::uint32_t device_id) { return impl::get_pci_bus_id(device_id); }

/**
 * @brief Get warp/subgroup size.
 * @note For SYCL implementation, this is fixed as a compile-time constant. However, developper
 * should profile kernel by picking one of the values in caribou::impl::get_subgroup_size. There is
 * no universal optimal value.
 */
inline std::uint32_t get_warp_size(void) { return impl::get_warp_size(); }

}  // namespace caribou
