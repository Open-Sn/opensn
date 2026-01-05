/*
 * Created on Fri, April 11 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <chrono>       // std::chrono::duration, std::milli
#include <memory>       // std::shared_ptr
#include <type_traits>  // std::remove_pointer_t

#include "api_mapping.hpp"  // GPU_API
#include "exception.hpp"    // caribou::check_error

namespace caribou {

// CUDA/HIP Event
// --------------

namespace impl {
using EventPtr = std::shared_ptr<std::remove_pointer<::GPU_API(Event_t)>::type>;
}  // namespace impl

struct InterprocedualEventHandle {
    ::GPU_API(IpcEventHandle_t) ptr;
};

class Event : public impl::EventPtr {
  public:
    /// @name Constructors
    /// {
    /**
     * @brief Default constructor.
     * @details Create an event with a given flag.
     * @param blocking_sync If true, the host thread calling the event synchronization is blocked until the event
     * occurs.
     * @param disable_timing If true, the event does not record timing data.
     * @param interprocess If true, the event can be converted to interprocess event. This flag must be enabled along
     * with disable timing.
     */
    inline Event(bool blocking_sync = false, bool disable_timing = true, bool interprocess = false) :
    impl::EventPtr(Event::create_event(blocking_sync, disable_timing, interprocess), ::GPU_API(EventDestroy)) {}
    /** @brief Take ownership from another pre-created CUDA/HIP event.*/
    inline Event(::GPU_API(Event_t) event_ptr) : impl::EventPtr(event_ptr, ::GPU_API(EventDestroy)) {}
    /** @brief Construct an event from an interprocess event handle.*/
    inline Event(const InterprocedualEventHandle & ipc_handle) :
    impl::EventPtr(Event::create_event_from_ipc(ipc_handle), ::GPU_API(EventDestroy)) {}

    /** @brief Query event status.*/
    inline bool is_completed() const {
        ::GPU_API(Error_t) status = ::GPU_API(EventQuery)(this->get());
        if (status == ::GPU_API(Success)) {
            return true;
        } else if (status == ::GPU_API(ErrorNotReady)) {
            return false;
        } else {
            check_error(status);
            return false;  // Unreachable
        }
    }
    /// @}

    /** @brief Type-cast operator to ``::cudaEvent_t`` or ``::hipEvent_t``.*/
    inline operator ::GPU_API(Event_t)(void) const { return this->get(); }

    /**
     * @brief Synchronize with an event.
     * @details Wait until all works recorded in the event are completed.
     *
     * If the event was created with the `blocking_sync` flag set to true, the host thread calling this function is
     * blocked until the event occurs. Otherwise, the host thread will be busy-waiting until the event has been
     * completed by the device.
     */
    inline void synchronize() const { check_error(::GPU_API(EventSynchronize)(this->get())); }

    /**
     * @brief Get interprocess event handle.
     * @details Get the interprocess event handle for this event. The event must be created with the `interprocess` and
     * `disable_timing` flags set to true.
     *
     * The generated handle can be copied to another process and used to open the event in that process.
     */
    inline InterprocedualEventHandle get_ipc_handle() const {
        ::GPU_API(IpcEventHandle_t) handle;
        check_error(::GPU_API(IpcGetEventHandle)(&handle, this->get()));
        return InterprocedualEventHandle{handle};
    }

  private:
    /**
     * @brief Create an event.
     * @param blocking_sync If true, the host thread calling the event synchronization is blocked until the event occurs.
     * @param disable_timing If true, the event does not record timing data.
     * @param interprocess If true, the event can be used for interprocess event.
     */
    static inline ::GPU_API(Event_t) create_event(bool blocking_sync, bool disable_timing, bool interprocess) {
        ::GPU_API(Event_t) result;
        unsigned int flag = 0;
        if (blocking_sync) {
            flag |= GPU_API(EventBlockingSync);
        }
        if (disable_timing) {
            flag |= GPU_API(EventDisableTiming);
        }
        if (interprocess) {
            flag |= GPU_API(EventInterprocess);
        }
        check_error(::GPU_API(EventCreateWithFlags)(&result, flag));
        return result;
    }

    /**
     * @brief Create an event from an interprocess event handle.
     * @param ipc_handle The interprocess event handle.
     */
    static inline ::GPU_API(Event_t) create_event_from_ipc(const InterprocedualEventHandle & ipc_handle) {
        ::GPU_API(Event_t) result;
        check_error(::GPU_API(IpcOpenEventHandle)(&result, ipc_handle.ptr));
        return result;
    }
};

/** @brief Measure the elapsed time between 2 events (resolution ~0.5ms).*/
inline std::chrono::duration<float, std::milli> operator-(const Event & end, const Event & start) {
    float time_in_ms = 0.0f;
    check_error(::GPU_API(EventElapsedTime)(&time_in_ms, start.get(), end.get()));
    return std::chrono::duration<float, std::milli>(time_in_ms);
}

}  // namespace caribou
