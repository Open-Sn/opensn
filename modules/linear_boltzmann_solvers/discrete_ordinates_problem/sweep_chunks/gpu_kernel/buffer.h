// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>

namespace opensn::gpu_kernel
{

/// Solver buffer for the sweep.
template <std::size_t ndofs>
class Buffer
{
public:
  /// Construct and zero-fill the buffer.
  constexpr Buffer()
  {
    A_.fill(0.0);
    b_.fill(0.0);
    s_.fill(0.0);
  }

  /// Get pointer to the linear system matrix.
  constexpr double* A() { return A_.data(); }

  /// Get poitner to the vector containing the angular flux.
  constexpr double* b() { return b_.data(); }

  /// Get pointer to the source vector.
  constexpr double* s() { return s_.data(); }

private:
  /// Linear system matrix.
  std::array<double, ndofs * ndofs> A_;
  /// Vector containing the angular flux.
  std::array<double, ndofs> b_;
  /// Source.
  std::array<double, ndofs> s_;
};

} // namespace opensn::gpu_kernel
