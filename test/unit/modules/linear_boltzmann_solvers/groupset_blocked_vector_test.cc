// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset_blocked_vector.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

using ScopedBlockOverride = GroupsetBlockedVector::ScopedBlockOverride;

/// Three blocks, the middle one empty, as on a groupset whose block a rank does not hold.
const std::vector<std::size_t> kBlockSizes = {3, 0, 5};

/// Fills `vector` with 1, 2, 3, ... in storage order.
void
FillWithIndices(GroupsetBlockedVector& vector)
{
  const auto data = vector.Data();
  for (std::size_t i = 0; i < data.size(); ++i)
    data[i] = 1.0 + static_cast<double>(i);
}

} // namespace

TEST(GroupsetBlockedVectorTest, ResizeSetsBlocksAndZeroes)
{
  GroupsetBlockedVector vector(kBlockSizes);
  ASSERT_EQ(vector.GetNumBlocks(), 3U);
  EXPECT_EQ(vector.Size(), 8U);
  EXPECT_FALSE(vector.Empty());
  EXPECT_EQ(vector.GetBlockOffset(0), 0U);
  EXPECT_EQ(vector.GetBlockOffset(1), 3U);
  EXPECT_EQ(vector.GetBlockOffset(2), 3U);
  for (std::size_t b = 0; b < kBlockSizes.size(); ++b)
  {
    EXPECT_EQ(vector.GetBlockSize(b), kBlockSizes[b]);
    EXPECT_EQ(vector.Block(b).size(), kBlockSizes[b]);
    EXPECT_EQ(vector.OwnedBlock(b).data(), vector.Data().data() + vector.GetBlockOffset(b));
  }
  for (const double value : std::as_const(vector).Data())
    EXPECT_EQ(value, 0.0);
  EXPECT_FALSE(vector.IsRedirected());
}

TEST(GroupsetBlockedVectorTest, GenerationChangesOnlyWithTheBlocks)
{
  GroupsetBlockedVector vector;
  EXPECT_EQ(vector.GetNumBlocks(), 0U);
  EXPECT_TRUE(vector.Empty());

  vector.Resize(kBlockSizes);
  const auto generation = vector.GetGeneration();
  const double* address = vector.Data().data();

  // The same block sizes keep the buffer and the generation, and zero the entries.
  FillWithIndices(vector);
  vector.Resize(kBlockSizes);
  EXPECT_EQ(vector.GetGeneration(), generation);
  EXPECT_EQ(vector.Data().data(), address);
  for (const double value : std::as_const(vector).Data())
    EXPECT_EQ(value, 0.0);

  // Moving entries between blocks changes the generation even at the same total size.
  vector.Resize({5, 0, 3});
  EXPECT_NE(vector.GetGeneration(), generation);
}

TEST(GroupsetBlockedVectorTest, WholeVectorOperationsKeepTheBuffer)
{
  GroupsetBlockedVector vector(kBlockSizes);
  const auto generation = vector.GetGeneration();
  const double* address = vector.Data().data();

  vector.Fill(2.0);
  vector.Scale(-1.5);
  for (const double value : std::as_const(vector).Data())
    EXPECT_EQ(value, -3.0);

  std::vector<double> values(vector.Size());
  for (std::size_t i = 0; i < values.size(); ++i)
    values[i] = 0.25 * static_cast<double>(i);
  vector.Assign(values);
  EXPECT_EQ(std::vector<double>(vector.Data().begin(), vector.Data().end()), values);

  // Assigning the vector its own buffer leaves it unchanged.
  vector.Assign(std::as_const(vector).Data());
  EXPECT_EQ(std::vector<double>(vector.Data().begin(), vector.Data().end()), values);

  EXPECT_EQ(vector.Data().data(), address);
  EXPECT_EQ(vector.GetGeneration(), generation);

  const std::vector<double> short_values(vector.Size() - 1, 0.0);
  EXPECT_THROW(vector.Assign(short_values), std::invalid_argument);
}

TEST(GroupsetBlockedVectorTest, WritableOverrideRedirectsTheBlock)
{
  GroupsetBlockedVector vector(kBlockSizes);
  FillWithIndices(vector);
  const std::vector<double> owned(vector.Data().begin(), vector.Data().end());
  std::vector<double> target(kBlockSizes[2], -1.0);
  {
    const ScopedBlockOverride guard(vector, 2, std::span<double>(target));
    EXPECT_TRUE(vector.IsRedirected());
    EXPECT_EQ(vector.Block(2).data(), target.data());
    EXPECT_EQ(std::as_const(vector).Block(2).data(), target.data());
    EXPECT_EQ(vector.OwnedBlock(2).data(), vector.OwnedBlock(0).data() + 3);
    // Other blocks are not redirected.
    EXPECT_EQ(vector.Block(0).data(), vector.OwnedBlock(0).data());

    // Writes through the block reach the target, not the owned storage.
    for (auto& value : vector.Block(2))
      value = 7.0;
    for (const double value : target)
      EXPECT_EQ(value, 7.0);

    // Whole-vector access would see stale owned storage, so it is rejected.
    EXPECT_THROW(vector.Data(), std::logic_error);
    EXPECT_THROW(std::as_const(vector).Data(), std::logic_error);
    EXPECT_THROW(vector.Fill(0.0), std::logic_error);
    EXPECT_THROW(vector.Scale(2.0), std::logic_error);
    EXPECT_THROW(vector.Assign(owned), std::logic_error);
    EXPECT_THROW(vector.Resize(kBlockSizes), std::logic_error);
  }
  EXPECT_FALSE(vector.IsRedirected());
  EXPECT_EQ(vector.Block(2).data(), vector.OwnedBlock(2).data());
  EXPECT_EQ(std::vector<double>(vector.Data().begin(), vector.Data().end()), owned);
}

TEST(GroupsetBlockedVectorTest, ReadOnlyOverrideRejectsMutableAccess)
{
  GroupsetBlockedVector vector(kBlockSizes);
  const std::vector<double> target(kBlockSizes[0], 4.0);
  {
    const ScopedBlockOverride guard(vector, 0, std::span<const double>(target));
    EXPECT_EQ(std::as_const(vector).Block(0).data(), target.data());
    EXPECT_THROW(vector.Block(0), std::logic_error);
    EXPECT_THROW(GroupsetBlockedSpan<double>{vector}, std::logic_error);
    const GroupsetBlockedSpan<const double> blocks(std::as_const(vector));
    EXPECT_EQ(blocks.Block(0).data(), target.data());
    // A read-only redirect of one block does not affect the others.
    EXPECT_EQ(vector.Block(2).data(), vector.OwnedBlock(2).data());
  }
  EXPECT_EQ(vector.Block(0).data(), vector.OwnedBlock(0).data());
}

TEST(GroupsetBlockedVectorTest, OverrideIsReleasedWhenAnExceptionUnwinds)
{
  GroupsetBlockedVector vector(kBlockSizes);
  std::vector<double> target(kBlockSizes[2], 0.0);
  EXPECT_THROW(
    {
      const ScopedBlockOverride guard(vector, 2, std::span<double>(target));
      throw std::runtime_error("failure inside an operator application");
    },
    std::runtime_error);
  EXPECT_FALSE(vector.IsRedirected());
  EXPECT_EQ(vector.Block(2).data(), vector.OwnedBlock(2).data());
  EXPECT_NO_THROW(vector.Data());
}

TEST(GroupsetBlockedVectorTest, NestedOverridesRestoreInReverseOrder)
{
  GroupsetBlockedVector vector(kBlockSizes);
  std::vector<double> outer(kBlockSizes[0], 1.0);
  const std::vector<double> inner(kBlockSizes[0], 2.0);
  {
    const ScopedBlockOverride outer_guard(vector, 0, std::span<double>(outer));
    {
      const ScopedBlockOverride inner_guard(vector, 0, std::span<const double>(inner));
      EXPECT_EQ(std::as_const(vector).Block(0).data(), inner.data());
      EXPECT_THROW(vector.Block(0), std::logic_error);
    }
    // The outer, writable redirect is back.
    EXPECT_TRUE(vector.IsRedirected());
    EXPECT_EQ(vector.Block(0).data(), outer.data());
  }
  EXPECT_FALSE(vector.IsRedirected());
  EXPECT_EQ(vector.Block(0).data(), vector.OwnedBlock(0).data());
}

TEST(GroupsetBlockedVectorTest, OverrideRejectsAMismatchedTarget)
{
  GroupsetBlockedVector vector(kBlockSizes);
  std::vector<double> wrong_size(kBlockSizes[2] + 1, 0.0);
  EXPECT_THROW(ScopedBlockOverride(vector, 2, std::span<double>(wrong_size)),
               std::invalid_argument);
  EXPECT_THROW(ScopedBlockOverride(vector, 3, std::span<double>(wrong_size)),
               std::invalid_argument);
  EXPECT_FALSE(vector.IsRedirected());
  EXPECT_EQ(vector.Block(2).data(), vector.OwnedBlock(2).data());
}

TEST(GroupsetBlockedVectorTest, EmptyBlocksCanBeRedirected)
{
  // A rank without cells has empty blocks; solvers still install their overrides there.
  GroupsetBlockedVector vector(std::vector<std::size_t>{0, 0});
  EXPECT_TRUE(vector.Empty());
  {
    const ScopedBlockOverride guard(vector, 1, std::span<double>());
    EXPECT_TRUE(vector.IsRedirected());
    EXPECT_TRUE(vector.Block(1).empty());
  }
  EXPECT_FALSE(vector.IsRedirected());
  EXPECT_TRUE(vector.Data().empty());
}

TEST(GroupsetBlockedSpanTest, BlocksOfABufferFollowTheLayout)
{
  // Three groupsets of unequal sizes: groups [0, 2), [2, 5), and [5, 7).
  constexpr std::uint64_t num_nodes = 4;
  constexpr unsigned int num_moments = 3;
  const FluxMomentLayout layout(num_nodes, num_moments, 7, {{0, 2}, {2, 3}, {5, 2}});
  std::vector<double> buffer(layout.GetNumLocalDOFs());
  for (std::size_t i = 0; i < buffer.size(); ++i)
    buffer[i] = static_cast<double>(i);

  const GroupsetBlockedSpan<double> blocks(buffer, layout);
  ASSERT_EQ(blocks.GetNumBlocks(), layout.GetNumBlocks());
  for (std::size_t b = 0; b < layout.GetNumBlocks(); ++b)
  {
    EXPECT_EQ(blocks.Block(b).data(), buffer.data() + layout.GetBlockOffset(b));
    EXPECT_EQ(blocks.Block(b).size(), layout.GetBlockSize(b));
  }

  // A block-relative index into a group's block reads the entry of the whole-vector index.
  const auto group_blocks = blocks.GetGroupBlocks(layout);
  ASSERT_EQ(group_blocks.size(), 7U);
  for (std::uint64_t n = 0; n < num_nodes; ++n)
    for (unsigned int m = 0; m < num_moments; ++m)
      for (unsigned int g = 0; g < 7; ++g)
        EXPECT_EQ(group_blocks[g][layout.MapBlockDOF(n, m, g)], buffer[layout.MapDOF(n, m, g)]);

  const GroupsetBlockedSpan<const double> read_only(blocks);
  EXPECT_EQ(read_only.Block(1).data(), blocks.Block(1).data());

  // Both views are copyable values: copies hold the same blocks.
  GroupsetBlockedSpan<double> copy;
  copy = blocks;
  GroupsetBlockedSpan<const double> read_only_copy;
  read_only_copy = read_only;
  EXPECT_EQ(copy.Block(2).data(), blocks.Block(2).data());
  EXPECT_EQ(read_only_copy.Block(2).data(), blocks.Block(2).data());

  std::vector<double> short_buffer(buffer.size() - 1);
  EXPECT_THROW(GroupsetBlockedSpan<double>(short_buffer, layout), std::invalid_argument);
  const FluxMomentLayout other_layout(num_nodes, num_moments, 7, {{0, 7}});
  EXPECT_THROW(blocks.GetGroupBlocks(other_layout), std::invalid_argument);
}

TEST(GroupsetBlockedSpanTest, CapturesTheCurrentBlocksOfAVector)
{
  GroupsetBlockedVector vector(kBlockSizes);
  std::vector<double> target(kBlockSizes[2], 0.0);
  const ScopedBlockOverride guard(vector, 2, std::span<double>(target));
  const GroupsetBlockedSpan<double> blocks(vector);
  ASSERT_EQ(blocks.GetNumBlocks(), kBlockSizes.size());
  EXPECT_EQ(blocks.Block(0).data(), vector.OwnedBlock(0).data());
  EXPECT_EQ(blocks.Block(2).data(), target.data());
}
