// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

namespace opensn
{

class FluxMomentLayout;
template <typename T>
class GroupsetBlockedSpan;

/**
 * Owner of a local flux-moment vector, such as the scalar flux or the source moments: one
 * contiguous buffer split into blocks, one per groupset, in the order of FluxMomentLayout.
 *
 * Block(b) returns the current storage of block b. That is the owned storage, unless a
 * ScopedBlockOverride redirects the block to other memory of the same size for the guard's
 * lifetime, for example to let a transport sweep write directly into a solver's vector. Code that
 * can run while a block is redirected must fetch the block with Block when it uses it and must not
 * keep the span across calls. A read-only override makes the mutable Block throw, so code that only
 * reads a block must use the const overload.
 *
 * Whole-vector access (Data, Fill, Scale, Assign) requires that no block is redirected and throws
 * std::logic_error otherwise, so it cannot silently read or write stale owned storage.
 *
 * Fill, Scale and Assign keep the buffer address. Only Resize can reallocate the buffer, and it
 * then increments GetGeneration, so holders of pointers into the owned storage can detect that they
 * are stale.
 *
 * Installing and releasing overrides is not thread safe: do it outside threaded regions.
 */
class GroupsetBlockedVector
{
public:
  class ScopedBlockOverride;

  GroupsetBlockedVector() = default;

  /// Allocates zero-filled blocks with the given entry counts.
  explicit GroupsetBlockedVector(const std::vector<std::size_t>& block_sizes);

  // Overrides point at the vector, so it can be neither copied nor moved.
  GroupsetBlockedVector(const GroupsetBlockedVector&) = delete;
  GroupsetBlockedVector& operator=(const GroupsetBlockedVector&) = delete;
  GroupsetBlockedVector(GroupsetBlockedVector&&) = delete;
  GroupsetBlockedVector& operator=(GroupsetBlockedVector&&) = delete;
  ~GroupsetBlockedVector() = default;

  /**
   * Sets the block entry counts and every entry to zero. The generation changes when the block
   * sizes or the buffer address change. Throws std::logic_error while a block is redirected.
   */
  void Resize(const std::vector<std::size_t>& block_sizes);

  std::size_t GetNumBlocks() const { return block_offset_.size() - 1; }

  /// Index in the whole vector of the first entry of block `block`.
  std::size_t GetBlockOffset(std::size_t block) const
  {
    assert(block < GetNumBlocks() && "GroupsetBlockedVector: block out of range.");
    return block_offset_[block];
  }

  std::size_t GetBlockSize(std::size_t block) const
  {
    assert(block < GetNumBlocks() && "GroupsetBlockedVector: block out of range.");
    return block_offset_[block + 1] - block_offset_[block];
  }

  /// Total entry count of all blocks.
  std::size_t Size() const { return data_.size(); }

  bool Empty() const { return data_.empty(); }

  /// Changes whenever Resize changes the block sizes or reallocates the buffer.
  std::uint64_t GetGeneration() const { return generation_; }

  /// True while any block is redirected.
  bool IsRedirected() const { return num_redirected_ > 0; }

  /// Current storage of block `block`. Throws std::logic_error under a read-only override.
  std::span<double> Block(std::size_t block)
  {
    assert(block < GetNumBlocks() && "GroupsetBlockedVector: block out of range.");
    const auto& redirect = redirect_[block];
    if (not redirect.active)
      return {data_.data() + block_offset_[block], GetBlockSize(block)};
    if (redirect.read_only)
      ThrowReadOnlyBlock(block);
    return {redirect.writable, GetBlockSize(block)};
  }

  /// Current storage of block `block`, for reading.
  std::span<const double> Block(std::size_t block) const
  {
    assert(block < GetNumBlocks() && "GroupsetBlockedVector: block out of range.");
    const auto& redirect = redirect_[block];
    if (not redirect.active)
      return {data_.data() + block_offset_[block], GetBlockSize(block)};
    return {redirect.readable, GetBlockSize(block)};
  }

  /// Owned storage of block `block`, whether or not the block is redirected.
  std::span<double> OwnedBlock(std::size_t block)
  {
    assert(block < GetNumBlocks() && "GroupsetBlockedVector: block out of range.");
    return {data_.data() + block_offset_[block], GetBlockSize(block)};
  }

  std::span<const double> OwnedBlock(std::size_t block) const
  {
    assert(block < GetNumBlocks() && "GroupsetBlockedVector: block out of range.");
    return {data_.data() + block_offset_[block], GetBlockSize(block)};
  }

  /// The whole owned buffer. Throws std::logic_error while a block is redirected.
  std::span<double> Data()
  {
    if (IsRedirected())
      ThrowRedirected("Data");
    return data_;
  }

  /// The whole owned buffer, for reading. Throws std::logic_error while a block is redirected.
  std::span<const double> Data() const
  {
    if (IsRedirected())
      ThrowRedirected("Data");
    return data_;
  }

  /**
   * The current blocks, for kernels that write them. Throws std::logic_error if a block is
   * redirected read-only. Valid until an override changes or Resize runs.
   */
  GroupsetBlockedSpan<double> Blocks();

  /// The current blocks, for kernels that only read them.
  GroupsetBlockedSpan<const double> ConstBlocks() const;

  /// Sets every entry. Throws std::logic_error while a block is redirected.
  void Fill(double value);

  /// Multiplies every entry by `factor`. Throws std::logic_error while a block is redirected.
  void Scale(double factor);

  /**
   * Copies `values` into the whole buffer; `values` may be the buffer itself, but must not
   * otherwise overlap it. Throws std::invalid_argument unless `values` has Size() entries, and
   * std::logic_error while a block is redirected.
   */
  void Assign(std::span<const double> values);

private:
  struct Redirect
  {
    bool active = false;
    /// Set explicitly: an empty target, as on a rank without cells, may have a null address.
    bool read_only = false;
    /// Null for a read-only override.
    double* writable = nullptr;
    const double* readable = nullptr;
  };

  [[noreturn]] static void ThrowReadOnlyBlock(std::size_t block);
  [[noreturn]] static void ThrowRedirected(const char* operation);

  std::vector<double> data_;
  /// Offset of each block, followed by the total entry count.
  std::vector<std::size_t> block_offset_{0};
  std::vector<Redirect> redirect_;
  std::size_t num_redirected_ = 0;
  std::uint64_t generation_ = 0;
};

/**
 * RAII guard that redirects one block of a GroupsetBlockedVector to external memory of the block's
 * size and restores the previous storage when it is destroyed, including during stack unwinding.
 * Guards on one block may nest; they must be released in reverse order, which scoping ensures.
 * The external memory must outlive the guard. Install and release guards outside threaded regions.
 */
class GroupsetBlockedVector::ScopedBlockOverride
{
public:
  /**
   * Redirects block `block` to the writable `target`. Throws std::invalid_argument unless
   * `target` has the block's size.
   */
  ScopedBlockOverride(GroupsetBlockedVector& vector, std::size_t block, std::span<double> target);

  /**
   * Redirects block `block` to the read-only `target`: the mutable Block throws while the guard is
   * active. Throws std::invalid_argument unless `target` has the block's size.
   */
  ScopedBlockOverride(GroupsetBlockedVector& vector,
                      std::size_t block,
                      std::span<const double> target);

  ScopedBlockOverride(const ScopedBlockOverride&) = delete;
  ScopedBlockOverride& operator=(const ScopedBlockOverride&) = delete;
  ScopedBlockOverride(ScopedBlockOverride&&) = delete;
  ScopedBlockOverride& operator=(ScopedBlockOverride&&) = delete;

  ~ScopedBlockOverride();

private:
  void Install(bool read_only, double* writable, const double* readable, std::size_t size);

  GroupsetBlockedVector& vector_;
  std::size_t block_;
  Redirect previous_;
};

/**
 * One span per groupset block of a flux-moment vector, in block order: a non-owning view that
 * kernels take instead of a whole vector, so they work on a GroupsetBlockedVector whose blocks may
 * be redirected as well as on any buffer in the flux-moment layout.
 *
 * Built from a GroupsetBlockedVector, it captures the current blocks, honoring overrides, so it
 * must not outlive the overrides active when it was built. `T` is `double` or `const double`.
 */
template <typename T>
class GroupsetBlockedSpan
{
  static_assert(std::is_same_v<std::remove_const_t<T>, double>,
                "GroupsetBlockedSpan holds double or const double.");

public:
  GroupsetBlockedSpan() = default;

  /// Current blocks of `vector`. Throws std::logic_error if a block is redirected read-only.
  explicit GroupsetBlockedSpan(GroupsetBlockedVector& vector)
    requires(not std::is_const_v<T>)
  {
    blocks_.reserve(vector.GetNumBlocks());
    for (std::size_t b = 0; b < vector.GetNumBlocks(); ++b)
      blocks_.push_back(vector.Block(b));
  }

  /// Current blocks of `vector`, for reading.
  explicit GroupsetBlockedSpan(const GroupsetBlockedVector& vector)
    requires(std::is_const_v<T>)
  {
    blocks_.reserve(vector.GetNumBlocks());
    for (std::size_t b = 0; b < vector.GetNumBlocks(); ++b)
      blocks_.push_back(vector.Block(b));
  }

  /**
   * Read-only view of a writable one. A template, so that it is not the copy constructor of
   * GroupsetBlockedSpan<double>, which stays implicitly copyable and movable.
   */
  template <typename U>
    requires(std::is_const_v<T> and std::is_same_v<U, double>)
  explicit GroupsetBlockedSpan(const GroupsetBlockedSpan<U>& other)
  {
    blocks_.reserve(other.GetNumBlocks());
    for (std::size_t b = 0; b < other.GetNumBlocks(); ++b)
      blocks_.push_back(other.Block(b));
  }

  /**
   * Blocks of `buffer`, a whole vector in the order of `layout`. Throws std::invalid_argument
   * unless `buffer` has layout.GetNumLocalDOFs() entries.
   */
  GroupsetBlockedSpan(std::span<T> buffer, const FluxMomentLayout& layout);

  std::size_t GetNumBlocks() const { return blocks_.size(); }

  std::span<T> Block(std::size_t block) const
  {
    assert(block < blocks_.size() && "GroupsetBlockedSpan: block out of range.");
    return blocks_[block];
  }

  /**
   * The block that holds each group of `layout`, indexed by group, for kernels that read several
   * groupsets. Throws std::invalid_argument unless the block count matches `layout`.
   */
  std::vector<std::span<T>> GetGroupBlocks(const FluxMomentLayout& layout) const;

private:
  std::vector<std::span<T>> blocks_;
};

extern template class GroupsetBlockedSpan<double>;
extern template class GroupsetBlockedSpan<const double>;

inline GroupsetBlockedSpan<double>
GroupsetBlockedVector::Blocks()
{
  return GroupsetBlockedSpan<double>(*this);
}

inline GroupsetBlockedSpan<const double>
GroupsetBlockedVector::ConstBlocks() const
{
  return GroupsetBlockedSpan<const double>(*this);
}

} // namespace opensn
