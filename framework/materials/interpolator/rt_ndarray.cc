// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/interpolator/rt_ndarray.h"

#include <algorithm>
#include <numeric>

namespace opensn
{

RT_NDArray::RT_NDArray(std::span<const std::uint32_t> shape, const std::size_t itemsize)
{
  Construct(shape, itemsize);
}

std::span<double>
RT_NDArray::At(std::span<const std::uint32_t> index)
{
  const auto item_offset = ComputeItemOffsetChecked(index);
  return std::span<double>{storage_.data() + item_offset * itemsize_, itemsize_};
}

std::span<const double>
RT_NDArray::At(std::span<const std::uint32_t> index) const
{
  const auto item_offset = ComputeItemOffsetChecked(index);
  return std::span<const double>{storage_.data() + item_offset * itemsize_, itemsize_};
}

void
RT_NDArray::Construct(std::span<const std::uint32_t> shape, const std::size_t itemsize)
{
  if (shape.size() > MAX_DIM)
    throw std::invalid_argument("Shape rank exceeds MAX_DIM.");
  if (itemsize == 0)
    throw std::invalid_argument("itemsize must be positive.");

  ndim_ = shape.size();
  itemsize_ = itemsize;
  shape_.fill(0);
  strides_.fill(0);

  std::copy(shape.begin(), shape.end(), shape_.begin());
  if (std::any_of(shape_.begin(),
                  shape_.begin() + static_cast<std::ptrdiff_t>(ndim_),
                  [](const std::uint32_t size) { return size == 0; }))
  {
    throw std::invalid_argument("Shape entries must be positive.");
  }

  std::size_t num_items = 1;
  for (std::uint32_t d = ndim_; d-- > 0;)
  {
    strides_[d] = num_items;
    num_items *= shape_[d];
  }

  storage_.assign(num_items * itemsize_, 0.0);
}

std::size_t
RT_NDArray::ComputeItemOffsetChecked(std::span<const std::uint32_t> index) const
{
  if (index.size() != ndim_)
  {
    throw std::invalid_argument("Number of indices " + std::to_string(index.size()) +
                                " not equal to shape rank " + std::to_string(ndim_) + ".");
  }

  std::size_t item_offset = 0;
  for (std::uint32_t d = 0; d < ndim_; ++d)
  {
    if (index[d] >= shape_[d])
    {
      throw std::out_of_range("Index " + std::to_string(d) + " out of range " +
                              std::to_string(index[d]) + " must be <" + std::to_string(shape_[d]) +
                              ".");
    }
    item_offset += index[d] * strides_[d];
  }

  return item_offset;
}

} // namespace opensn
