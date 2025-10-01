/*
 * Created on Fri, April 11 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <chrono>       // std::chrono::duration, std::milli
#include <memory>       // std::shared_ptr
#include <type_traits>  // std::remove_pointer_t

#include "exception.hpp"  // cuda::check_cuda_error

namespace caribou {

// CUDA Event
// ----------

namespace cuda {
using EventImpl = std::shared_ptr<std::remove_pointer<::cudaEvent_t>::type>;
}  // namespace cuda

struct InterprocedualEventHandle {
    ::cudaIpcEventHandle_t ptr;
};

class Event : public cuda::EventImpl {
  public:
    /**
     * @brief Default constructor.
     * @details Create an event with a given flag.
     * @param blocking_sync If true, the host thread calling the event synchronization is blocked until the event
     * occurs.
     * @param disable_timing If true, the event does not record timing data.
     * @param interprocess If true, the event can be converted to interprocess event. This flag must be enabled along with
     * disable timing.
     */
    inline Event(bool blocking_sync = false, bool disable_timing = true, bool interprocess = false) :
    cuda::EventImpl(Event::create_event(blocking_sync, disable_timing, interprocess), ::cudaEventDestroy) {}
    /** @brief Take ownership from another pre-created CUDA event.*/
    inline Event(::cudaEvent_t event_ptr) : cuda::EventImpl(event_ptr, ::cudaEventDestroy) {}
    /** @brief Construct an event from an interprocess event handle.*/
    inline Event(const InterprocedualEventHandle & ipc_handle) :
    cuda::EventImpl(Event::create_event_from_ipc(ipc_handle), ::cudaEventDestroy) {}

    /** @brief Query event status.*/
    inline bool is_completed() const {
        ::cudaError_t status = ::cudaEventQuery(this->get());
        if (status == ::cudaSuccess) {
            return true;
        } else if (status == ::cudaErrorNotReady) {
            return false;
        } else {
            cuda::check_cuda_error(status);
            return false;  // Unreachable
        }
    }

    /**
     * @brief Synchronize with an event.
     * @details Wait until all works recorded in the event are completed.
     *
     * If the event was created with the `blocking_sync` flag set to true, the host thread calling this function is
     * blocked until the event occurs. Otherwise, the host thread will be busy-waiting until the event has been
     * completed by the device.
     */
    inline void synchronize() const { cuda::check_cuda_error(::cudaEventSynchronize(this->get())); }

    /**
     * @brief Get interprocess event handle.
     * @details Get the interprocess event handle for this event. The event must be created with the `interprocess` and
     * `disable_timing` flags set to true.
     *
     * The generated handle can be copied to another process and used to open the event in that process.
     */
    inline InterprocedualEventHandle get_ipc_handle() const {
        ::cudaIpcEventHandle_t handle;
        cuda::check_cuda_error(::cudaIpcGetEventHandle(&handle, this->get()));
        return InterprocedualEventHandle{handle};
    }

  private:
    /**
     * @brief Create an event.
     * @param blocking_sync If true, the host thread calling the event synchronization is blocked until the event occurs.
     * @param disable_timing If true, the event does not record timing data.
     * @param interprocess If true, the event can be used for interprocess event.
     */
    static inline ::cudaEvent_t create_event(bool blocking_sync, bool disable_timing, bool interprocess) {
        ::cudaEvent_t result;
        unsigned int flag = 0;
        if (blocking_sync) {
            flag |= cudaEventBlockingSync;
        }
        if (disable_timing) {
            flag |= cudaEventDisableTiming;
        }
        if (interprocess) {
            flag |= cudaEventInterprocess;
        }
        cuda::check_cuda_error(::cudaEventCreateWithFlags(&result, flag));
        return result;
    }

    /**
     * @brief Create an event from an interprocess event handle.
     * @param ipc_handle The interprocess event handle.
     */
    static inline ::cudaEvent_t create_event_from_ipc(const InterprocedualEventHandle & ipc_handle) {
        ::cudaEvent_t result;
        cuda::check_cuda_error(::cudaIpcOpenEventHandle(&result, ipc_handle.ptr));
        return result;
    }
};

/** @brief Measure the elapsed time between 2 events (resolution ~0.5ms).*/
inline std::chrono::duration<float, std::milli> operator-(const Event & end, const Event & start) {
    float time_in_ms = 0.0f;
    cuda::check_cuda_error(::cudaEventElapsedTime(&time_in_ms, start.get(), end.get()));
    return std::chrono::duration<float, std::milli>(time_in_ms);
}

}  // namespace caribou
