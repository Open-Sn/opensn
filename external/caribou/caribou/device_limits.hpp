// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <optional>

namespace caribou {

/** @brief Device-wide limits for one-dimensional blocks. */
struct LaunchLimits {
    /** @brief Maximum threads in a workgroup. */
    std::size_t max_block_threads;
    /** @brief Maximum threads along logical block x (SYCL local dimension 2). */
    std::size_t max_block_x;
    /** @brief Parameter capacity in bytes, or no value when the backend has no query. */
    std::optional<std::size_t> max_parameter_bytes;
};

}  // namespace caribou
