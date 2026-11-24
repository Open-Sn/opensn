// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/gpu/device_binding.h"
#include <exception>
#include <thread>
#include <vector>

namespace opensn
{

/// Run `function(i)` for i in [0, count) across `num_threads` threads, strided.
/// Exceptions thrown by any worker are propagated to the caller (first wins).
template <typename Function>
void
ParallelFor(size_t count, size_t num_threads, Function function)
{
  std::vector<std::exception_ptr> exceptions(num_threads);
  {
    std::vector<std::jthread> workers;
    workers.reserve(num_threads);
    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id)
      workers.emplace_back(
        [&, thread_id]()
        {
          try
          {
            UseBoundDevice();
            for (size_t i = thread_id; i < count; i += num_threads)
              function(i);
          }
          catch (...)
          {
            exceptions[thread_id] = std::current_exception();
          }
        });
  }

  for (const auto& exception : exceptions)
    if (exception)
      std::rethrow_exception(exception);
}

} // namespace opensn
