// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/device/dof_limits.h"
#include "caribou/main.hpp"
#include <array>
#include <cstddef>

namespace opensn::gpu_kernel
{

template <std::size_t ndofs>
class MatBuffer;

/// Solver buffer for the sweep.
template <std::size_t ndofs>
  requires(ndofs <= max_dof_gpu_register)
class MatBuffer<ndofs>
{
public:
  constexpr MatBuffer(double* un_used) { A_.fill(0.0); }

  struct row_pointer
  {
  public:
    constexpr row_pointer(double* data) : data_(data) {}
    constexpr double& operator[](std::uint32_t idx) { return data_[idx]; }
    constexpr const double& operator[](std::uint32_t idx) const { return data_[idx]; }
    constexpr row_pointer operator+(std::uint32_t offset) const
    {
      return row_pointer(data_ + ndofs * offset);
    }
    constexpr row_pointer operator++()
    {
      data_ += ndofs;
      return *this;
    }
    constexpr row_pointer operator++(int)
    {
      row_pointer temp = *this;
      data_ += ndofs;
      return temp;
    }

  private:
    double* data_;
  };
  /// Get row pointer for the linear system matrix.
  constexpr row_pointer row(std::uint32_t idx) { return row_pointer(A_.data() + idx * ndofs); }

  struct element_pointer
  {
  public:
    constexpr element_pointer(double* data) : data_(data) {}
    constexpr double& operator*() { return *data_; }
    constexpr const double& operator*() const { return *data_; }
    constexpr element_pointer operator++()
    {
      ++data_;
      return *this;
    }
    constexpr element_pointer operator++(int)
    {
      element_pointer temp = *this;
      ++data_;
      return temp;
    }

  private:
    double* data_;
  };
  /// Get element pointer for the linear system matrix.
  constexpr element_pointer element(std::uint32_t row_idx = 0, std::uint32_t col_idx = 0)
  {
    return element_pointer(A_.data() + row_idx * ndofs + col_idx);
  }

private:
  /// Linear system matrix.
  std::array<double, ndofs * ndofs> A_;
};

template <std::size_t ndofs>
  requires(ndofs > max_dof_gpu_register)
class MatBuffer<ndofs>
{
public:
  __device__ MatBuffer(double* data) : A_(data)
  {
    _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs * ndofs; ++i)
    {
      *data = 0.0;
      data += warpSize;
    }
  }

  struct row_pointer
  {
  public:
    constexpr row_pointer(double* data) : data_(data) {}
    constexpr double& operator[](std::uint32_t idx) { return data_[idx * warpSize]; }
    constexpr const double& operator[](std::uint32_t idx) const { return data_[idx * warpSize]; }
    constexpr row_pointer operator+(std::uint32_t offset) const
    {
      return row_pointer(data_ + ndofs * warpSize * offset);
    }
    constexpr row_pointer operator++()
    {
      data_ += ndofs * warpSize;
      return *this;
    }
    constexpr row_pointer operator++(int)
    {
      row_pointer temp = *this;
      data_ += ndofs * warpSize;
      return temp;
    }

  private:
    double* data_;
  };
  /// Get row pointer for the linear system matrix.
  constexpr row_pointer row(std::uint32_t idx) { return row_pointer(A_ + idx * ndofs * warpSize); }

  struct element_pointer
  {
  public:
    constexpr element_pointer(double* data) : data_(data) {}
    constexpr double& operator*() { return *data_; }
    constexpr const double& operator*() const { return *data_; }
    constexpr element_pointer operator++()
    {
      data_ += warpSize;
      return *this;
    }
    constexpr element_pointer operator++(int)
    {
      element_pointer temp = *this;
      data_ += warpSize;
      return temp;
    }

  private:
    double* data_;
  };
  /// Get element pointer for the linear system matrix.
  constexpr element_pointer element(std::uint32_t row_idx = 0, std::uint32_t col_idx = 0)
  {
    return element_pointer(A_ + (row_idx * ndofs + col_idx) * warpSize);
  }

private:
  double* A_ = nullptr;
};

} // namespace opensn::gpu_kernel
