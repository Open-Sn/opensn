// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset_blocked_vector.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <string>
#include <utility>

namespace opensn
{

GroupsetBlockedVector::GroupsetBlockedVector(const std::vector<std::size_t>& block_sizes)
{
  Resize(block_sizes);
}

void
GroupsetBlockedVector::Resize(const std::vector<std::size_t>& block_sizes)
{
  if (IsRedirected())
    ThrowRedirected("Resize");

  std::vector<std::size_t> block_offset(1, 0);
  block_offset.reserve(block_sizes.size() + 1);
  for (const auto size : block_sizes)
    block_offset.push_back(block_offset.back() + size);

  const double* old_address = data_.data();
  data_.assign(block_offset.back() + GetTailCapacity(), 0.0);
  if (block_offset != block_offset_ or data_.data() != old_address)
    ++generation_;
  block_offset_ = std::move(block_offset);
  redirect_.assign(block_sizes.size(), Redirect{});
}

void
GroupsetBlockedVector::ReserveTail(std::size_t num_entries)
{
  if (IsRedirected())
    ThrowRedirected("ReserveTail");
  if (num_entries <= GetTailCapacity())
    return;
  const double* old_address = data_.data();
  // Reallocate explicitly rather than through resize, whose growth policy may reserve more.
  std::vector<double> data(Size() + num_entries, 0.0);
  std::copy_n(data_.begin(), Size(), data.begin());
  data_ = std::move(data);
  if (data_.data() != old_address)
    ++generation_;
}

std::span<double>
GroupsetBlockedVector::OwnedBlockWithTail(std::size_t block, std::size_t tail_size)
{
  OpenSnInvalidArgumentIf(block >= GetNumBlocks(),
                          "GroupsetBlockedVector::OwnedBlockWithTail: block " +
                            std::to_string(block) + " does not exist; the vector has " +
                            std::to_string(GetNumBlocks()) + " blocks.");
  OpenSnInvalidArgumentIf(tail_size > 0 and block + 1 != GetNumBlocks(),
                          "GroupsetBlockedVector::OwnedBlockWithTail: only the last block can "
                          "extend into the tail, but block " +
                            std::to_string(block) + " of " + std::to_string(GetNumBlocks()) +
                            " was given a tail of " + std::to_string(tail_size) + " entries.");
  OpenSnInvalidArgumentIf(tail_size > GetTailCapacity(),
                          "GroupsetBlockedVector::OwnedBlockWithTail: requested " +
                            std::to_string(tail_size) + " tail entries, but the tail has " +
                            std::to_string(GetTailCapacity()) + ".");
  return {data_.data() + block_offset_[block], GetBlockSize(block) + tail_size};
}

void
GroupsetBlockedVector::Fill(double value)
{
  const auto data = Data();
  std::fill(data.begin(), data.end(), value);
}

void
GroupsetBlockedVector::Scale(double factor)
{
  for (auto& value : Data())
    value *= factor;
}

void
GroupsetBlockedVector::Assign(std::span<const double> values)
{
  const auto data = Data();
  OpenSnInvalidArgumentIf(values.size() != data.size(),
                          "GroupsetBlockedVector::Assign: got " + std::to_string(values.size()) +
                            " values, expected " + std::to_string(data.size()) + ".");
  // Assigning the vector its own buffer leaves it unchanged.
  if (values.data() != data.data())
    std::copy(values.begin(), values.end(), data.begin());
}

void
GroupsetBlockedVector::ThrowReadOnlyBlock(std::size_t block)
{
  OpenSnLogicalError("GroupsetBlockedVector: block " + std::to_string(block) +
                     " is redirected read-only and cannot be written.");
}

void
GroupsetBlockedVector::ThrowRedirected(const char* operation)
{
  OpenSnLogicalError(std::string("GroupsetBlockedVector::") + operation +
                     " needs the whole owned vector, but a block is redirected.");
}

GroupsetBlockedVector::ScopedBlockOverride::ScopedBlockOverride(GroupsetBlockedVector& vector,
                                                                std::size_t block,
                                                                std::span<double> target)
  : vector_(vector), block_(block)
{
  Install(false, target.data(), target.data(), target.size());
}

GroupsetBlockedVector::ScopedBlockOverride::ScopedBlockOverride(GroupsetBlockedVector& vector,
                                                                std::size_t block,
                                                                std::span<const double> target)
  : vector_(vector), block_(block)
{
  Install(true, nullptr, target.data(), target.size());
}

void
GroupsetBlockedVector::ScopedBlockOverride::Install(bool read_only,
                                                    double* writable,
                                                    const double* readable,
                                                    std::size_t size)
{
  OpenSnInvalidArgumentIf(block_ >= vector_.GetNumBlocks(),
                          "ScopedBlockOverride: block " + std::to_string(block_) +
                            " does not exist; the vector has " +
                            std::to_string(vector_.GetNumBlocks()) + " blocks.");
  OpenSnInvalidArgumentIf(size != vector_.GetBlockSize(block_),
                          "ScopedBlockOverride: block " + std::to_string(block_) + " has " +
                            std::to_string(vector_.GetBlockSize(block_)) +
                            " entries, but the target has " + std::to_string(size) + ".");
  auto& redirect = vector_.redirect_[block_];
  previous_ = redirect;
  if (not previous_.active)
    ++vector_.num_redirected_;
  redirect = {true, read_only, writable, readable};
}

GroupsetBlockedVector::ScopedBlockOverride::~ScopedBlockOverride()
{
  vector_.redirect_[block_] = previous_;
  if (not previous_.active)
    --vector_.num_redirected_;
}

template <typename T>
GroupsetBlockedSpan<T>::GroupsetBlockedSpan(std::span<T> buffer, const FluxMomentLayout& layout)
{
  OpenSnInvalidArgumentIf(buffer.size() != layout.GetNumLocalDOFs(),
                          "GroupsetBlockedSpan: the buffer has " + std::to_string(buffer.size()) +
                            " entries, but the flux-moment layout has " +
                            std::to_string(layout.GetNumLocalDOFs()) + ".");
  blocks_.reserve(layout.GetNumBlocks());
  for (std::size_t b = 0; b < layout.GetNumBlocks(); ++b)
    blocks_.push_back(buffer.subspan(layout.GetBlockOffset(b), layout.GetBlockSize(b)));
}

template <typename T>
std::vector<std::span<T>>
GroupsetBlockedSpan<T>::GetGroupBlocks(const FluxMomentLayout& layout) const
{
  OpenSnInvalidArgumentIf(blocks_.size() != layout.GetNumBlocks(),
                          "GroupsetBlockedSpan: has " + std::to_string(blocks_.size()) +
                            " blocks, but the flux-moment layout has " +
                            std::to_string(layout.GetNumBlocks()) + ".");
  std::vector<std::span<T>> group_blocks(layout.GetNumGroups());
  for (unsigned int g = 0; g < layout.GetNumGroups(); ++g)
    group_blocks[g] = blocks_[layout.GetGroupBlock(g)];
  return group_blocks;
}

template class GroupsetBlockedSpan<double>;
template class GroupsetBlockedSpan<const double>;

} // namespace opensn
