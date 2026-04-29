/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#include "caribou/sycl/device.hpp"  // caribou::impl::gpus

namespace caribou::impl {

std::vector<::sycl::device> gpus(::sycl::device::get_devices(::sycl::info::device_type::gpu));

}  // namespace caribou::impl
