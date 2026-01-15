// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>

namespace opensn::cbc_gpu_kernel
{

/// Linear system buffer for the device CBC sweep chunk kernel.
template <std::size_t ndofs>
class LinearSystemBuffer
{
public:
  /// Constructor and zero-fill the buffer.
  constexpr LinearSystemBuffer()
  {
    A_.fill(0.0);
    b_.fill(0.0);
    s_.fill(0.0);
  }

  /// Get pointer to linear system matrix.
  constexpr double* A() { return A_.data(); }

  /// Get pointer to angular flux vector.
  constexpr double* b() { return b_.data(); }

  /// Get pointer to source vector.
  constexpr double* s() { return s_.data(); }

private:
  /// Linear system matrix.
  std::array<double, ndofs * ndofs> A_;

  /// Angular flux vector.
  std::array<double, ndofs> b_;

  /// Source vector.
  std::array<double, ndofs> s_;
};

} // namespace opensn::cbc_gpu_kernel