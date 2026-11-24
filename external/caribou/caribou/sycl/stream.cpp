/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#include "caribou/sycl/stream.hpp"

#include <atomic>    // std::atomic, std::memory_order_acquire, std::memory_order_release
#include <iostream>  // std::clog
#include <limits>    // std::numeric_limits
#include <memory>    // std::unique_ptr
#include <mutex>     // std::call_once, std::once_flag
#include <optional>  // std::optional
#include <sstream>   // std::ostringstream

#include "caribou/sycl/device.hpp"  // caribou::impl::get_current_device, caribou::impl::get_gpu_devices,
                                    // caribou::impl::get_num_gpus

namespace caribou::impl {

namespace {

/** @brief Default (null) queue of one GPU, created when the GPU first requests it.*/
struct DefaultQueue {
    std::once_flag created;
    std::atomic<bool> ready{false};
    std::optional<StreamNativeHandle> queue;
};

/** @brief Default queues of all GPUs, indexed like caribou::impl::get_gpu_devices.*/
DefaultQueue * get_default_queues(void) {
    static const std::unique_ptr<DefaultQueue[]> queues(new DefaultQueue[get_num_gpus()]);
    return queues.get();
}

}  // namespace

StreamNativeHandle & get_null_stream_handle(void) {
    if (get_num_gpus() == 0) {
        throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::runtime), "No GPU device found.");
    }
    // Only the requested GPU must support fp64, so its queue is created and checked on its own.
    const std::uint32_t device_id = get_current_device();
    DefaultQueue & entry = get_default_queues()[device_id];
    std::call_once(entry.created, [&entry, device_id]() {
        const ::sycl::device & gpu_device = get_gpu_devices()[device_id];
        if (not gpu_device.has(::sycl::aspect::fp64)) {
            std::ostringstream message;
            message << "GPU device " << device_id << " (\""
                    << gpu_device.get_info<::sycl::info::device::name>() << "\") doesn't support fp64.";
            throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::runtime), message.str());
        }
        entry.queue.emplace(gpu_device, ::sycl::property::queue::in_order{});
        entry.ready.store(true, std::memory_order_release);
    });
    return *entry.queue;
}

StreamNativeHandle & get_owner_stream_handle(const void * ptr) {
    // Every allocation is made on a default queue. The current GPU usually made it, so its queue is
    // tried first.
    const std::uint32_t num_gpus = get_num_gpus();
    const std::uint32_t current_device = get_current_device();
    DefaultQueue * queues = get_default_queues();
    for (std::uint32_t offset = 0; offset < num_gpus; ++offset) {
        DefaultQueue & entry = queues[(current_device + offset) % num_gpus];
        if (entry.ready.load(std::memory_order_acquire) &&
            ::sycl::get_pointer_type(ptr, entry.queue->get_context()) != ::sycl::usm::alloc::unknown) {
            return *entry.queue;
        }
    }
    throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::invalid),
                            "Memory was not allocated by caribou on any GPU.");
}

std::uint32_t get_stream_device(const Stream & stream) {
    const std::vector<::sycl::device> & gpus = get_gpu_devices();
    const ::sycl::device stream_device = stream.get_device();
    for (std::uint32_t index = 0; index < gpus.size(); ++index) {
        if (gpus[index] == stream_device) {
            return index;
        }
    }
    std::clog << "Warning: stream is not associated with a GPU (it targets \""
              << stream_device.get_info<::sycl::info::device::name>() << "\").\n";
    return std::numeric_limits<std::uint32_t>::max();
}

}  // namespace caribou::impl
