// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/flux_moment_layout.h"
#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{

constexpr std::uint64_t kNodes = 4;
constexpr unsigned int kMoments = 3;
constexpr unsigned int kGroups = 7;

/// Three groupsets of unequal sizes: groups [0, 2), [2, 5), and [5, 7).
const std::vector<FluxMomentLayout::GroupRange> kGroupsets = {{0, 2}, {2, 3}, {5, 2}};

FluxMomentLayout
MakeLayout()
{
  return {kNodes, kMoments, kGroups, kGroupsets};
}

} // namespace

TEST(FluxMomentLayoutTest, MapsEveryEntryExactlyOnce)
{
  const auto layout = MakeLayout();
  ASSERT_EQ(layout.GetNumLocalDOFs(), kNodes * kMoments * kGroups);

  std::vector<int> hits(layout.GetNumLocalDOFs(), 0);
  for (std::uint64_t n = 0; n < kNodes; ++n)
    for (unsigned int m = 0; m < kMoments; ++m)
      for (unsigned int g = 0; g < kGroups; ++g)
      {
        const auto index = layout.MapDOF(n, m, g);
        ASSERT_LT(index, hits.size());
        ++hits[index];
      }
  for (const int count : hits)
    EXPECT_EQ(count, 1);
}

TEST(FluxMomentLayoutTest, GroupsOfAGroupsetAreContiguous)
{
  const auto layout = MakeLayout();
  for (const auto& range : kGroupsets)
    for (std::uint64_t n = 0; n < kNodes; ++n)
      for (unsigned int m = 0; m < kMoments; ++m)
      {
        const auto first = layout.MapDOF(n, m, range.first_group);
        for (unsigned int k = 0; k < range.num_groups; ++k)
          EXPECT_EQ(layout.MapDOF(n, m, range.first_group + k), first + k);
      }
}

TEST(FluxMomentLayoutTest, MomentAndNodeStrides)
{
  const auto layout = MakeLayout();
  for (unsigned int g = 0; g < kGroups; ++g)
  {
    const auto stride = layout.GetMomentStride(g);
    for (std::uint64_t n = 0; n < kNodes; ++n)
      for (unsigned int m = 0; m + 1 < kMoments; ++m)
        EXPECT_EQ(layout.MapDOF(n, m + 1, g) - layout.MapDOF(n, m, g), stride);
    for (std::uint64_t n = 0; n + 1 < kNodes; ++n)
      EXPECT_EQ(layout.MapDOF(n + 1, 0, g) - layout.MapDOF(n, 0, g), kMoments * stride);
  }
}

TEST(FluxMomentLayoutTest, OrderIsGroupsetNodeMomentGroup)
{
  const auto layout = MakeLayout();
  EXPECT_FALSE(layout.IsNodeMomentGroupOrder());
  std::size_t block_offset = 0;
  for (const auto& range : kGroupsets)
  {
    for (std::uint64_t n = 0; n < kNodes; ++n)
      for (unsigned int m = 0; m < kMoments; ++m)
        for (unsigned int k = 0; k < range.num_groups; ++k)
          EXPECT_EQ(layout.MapDOF(n, m, range.first_group + k),
                    block_offset + (n * kMoments + m) * range.num_groups + k);
    block_offset += kNodes * kMoments * range.num_groups;
  }
}

TEST(FluxMomentLayoutTest, BlocksFollowTheGroupsets)
{
  const auto layout = MakeLayout();
  ASSERT_EQ(layout.GetNumBlocks(), kGroupsets.size());
  std::size_t block_offset = 0;
  for (std::size_t b = 0; b < kGroupsets.size(); ++b)
  {
    const auto& range = kGroupsets[b];
    EXPECT_EQ(layout.GetBlockOffset(b), block_offset);
    EXPECT_EQ(layout.GetBlockSize(b), kNodes * kMoments * range.num_groups);
    EXPECT_EQ(layout.GetBlockSizes()[b], layout.GetBlockSize(b));
    for (unsigned int k = 0; k < range.num_groups; ++k)
      EXPECT_EQ(layout.GetGroupBlock(range.first_group + k), b);
    block_offset += layout.GetBlockSize(b);
  }
  EXPECT_EQ(block_offset, layout.GetNumLocalDOFs());
}

TEST(FluxMomentLayoutTest, BlockIndexIsRelativeToTheGroupsetBlock)
{
  const auto layout = MakeLayout();
  for (const auto& range : kGroupsets)
    for (std::uint64_t n = 0; n < kNodes; ++n)
      for (unsigned int m = 0; m < kMoments; ++m)
        for (unsigned int k = 0; k < range.num_groups; ++k)
        {
          const auto g = range.first_group + k;
          const auto block_index = layout.MapBlockDOF(n, m, g);
          EXPECT_EQ(block_index, (n * kMoments + m) * range.num_groups + k);
          EXPECT_EQ(layout.GetBlockOffset(layout.GetGroupBlock(g)) + block_index,
                    layout.MapDOF(n, m, g));
        }
}

TEST(FluxMomentLayoutTest, SingleGroupsetIsNodeMomentGroup)
{
  const FluxMomentLayout layout(kNodes, kMoments, kGroups, {{0, kGroups}});
  EXPECT_TRUE(layout.IsNodeMomentGroupOrder());
  for (std::uint64_t n = 0; n < kNodes; ++n)
    for (unsigned int m = 0; m < kMoments; ++m)
      for (unsigned int g = 0; g < kGroups; ++g)
        EXPECT_EQ(layout.MapDOF(n, m, g), (n * kMoments + m) * kGroups + g);
}

TEST(FluxMomentLayoutTest, NodeMomentGroupOrderRoundTrip)
{
  const auto layout = MakeLayout();
  std::vector<double> ordered(layout.GetNumLocalDOFs());
  for (std::size_t k = 0; k < ordered.size(); ++k)
    ordered[k] = 0.5 + static_cast<double>(k);

  const auto mapped = layout.FromNodeMomentGroupOrder(ordered);
  for (std::uint64_t n = 0; n < kNodes; ++n)
    for (unsigned int m = 0; m < kMoments; ++m)
      for (unsigned int g = 0; g < kGroups; ++g)
        EXPECT_EQ(mapped[layout.MapDOF(n, m, g)], ordered[(n * kMoments + m) * kGroups + g]);

  EXPECT_EQ(layout.ToNodeMomentGroupOrder(mapped), ordered);
}

TEST(FluxMomentLayoutTest, OrderConversionRejectsWrongSize)
{
  const auto layout = MakeLayout();
  const std::vector<double> short_values(layout.GetNumLocalDOFs() - 1, 0.0);
  EXPECT_THROW(layout.ToNodeMomentGroupOrder(short_values), std::invalid_argument);
  EXPECT_THROW(layout.FromNodeMomentGroupOrder(short_values), std::invalid_argument);
}

TEST(FluxMomentLayoutTest, ZeroLocalNodes)
{
  const FluxMomentLayout layout(0, kMoments, kGroups, kGroupsets);
  EXPECT_EQ(layout.GetNumLocalDOFs(), 0U);
  ASSERT_EQ(layout.GetNumBlocks(), kGroupsets.size());
  for (std::size_t b = 0; b < kGroupsets.size(); ++b)
  {
    EXPECT_EQ(layout.GetBlockOffset(b), 0U);
    EXPECT_EQ(layout.GetBlockSize(b), 0U);
  }
  EXPECT_TRUE(layout.ToNodeMomentGroupOrder({}).empty());
  EXPECT_TRUE(layout.FromNodeMomentGroupOrder({}).empty());
}

TEST(FluxMomentLayoutTest, DefaultHasNoBlocks)
{
  const FluxMomentLayout layout;
  EXPECT_EQ(layout.GetNumBlocks(), 0U);
  EXPECT_TRUE(layout.GetBlockSizes().empty());
}

TEST(FluxMomentLayoutTest, RejectsGroupsetsThatDoNotPartitionTheGroups)
{
  using Ranges = std::vector<FluxMomentLayout::GroupRange>;
  // No groupset, a first groupset after group 0, a gap, an overlap, an empty groupset, a groupset
  // past the last group, and groupsets that stop short of the last group.
  const std::vector<Ranges> invalid = {
    {}, {{1, 6}}, {{0, 2}, {3, 4}}, {{0, 3}, {2, 5}}, {{0, 0}, {0, 7}}, {{0, 8}}, {{0, 2}, {2, 3}}};
  for (const auto& ranges : invalid)
    EXPECT_THROW(FluxMomentLayout(kNodes, kMoments, kGroups, ranges), std::logic_error);
}
