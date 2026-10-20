// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_message_struct.h"
#include "framework/mpi/sweep_communicator.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{

/// Splits the messages of one location with `num_unknowns` doubles and checks that the blocks
/// partition [0, num_unknowns) into nonempty, contiguous blocks inside the buffer, each of at
/// most `max_message_bytes` bytes.
void
ExpectValidSplit(const SweepCommunicator& sweep_communicator,
                 const std::size_t num_unknowns,
                 const int max_message_bytes)
{
  const std::vector<int> locations = {mpi_comm.rank()};
  std::vector<std::vector<AAH_MessageDetails>> msg_data;
  std::vector<std::vector<bool>> msg_received;
  int max_num_messages = 0;
  SetupMessageData(
    locations,
    [num_unknowns](std::size_t) { return num_unknowns; },
    msg_data,
    &msg_received,
    max_num_messages,
    sweep_communicator,
    max_message_bytes);

  ASSERT_EQ(msg_data.size(), 1U);
  ASSERT_EQ(msg_received.size(), 1U);
  const auto& blocks = msg_data.front();
  ASSERT_FALSE(blocks.empty()) << "num_unknowns=" << num_unknowns;

  long long next_position = 0;
  long long min_size = blocks.front().size;
  long long max_bytes = 0;
  bool contiguous = true;
  bool inside_buffer = true;
  for (const auto& block : blocks)
  {
    EXPECT_EQ(block.peer, mpi_comm.rank());
    contiguous = contiguous and (block.block_pos == next_position);
    inside_buffer =
      inside_buffer and block.block_pos >= 0 and
      static_cast<long long>(block.block_pos) + block.size <= static_cast<long long>(num_unknowns);
    next_position += block.size;
    min_size = std::min<long long>(min_size, block.size);
    max_bytes = std::max(
      max_bytes, static_cast<long long>(block.size) * static_cast<long long>(sizeof(double)));
  }

  EXPECT_GE(min_size, 1) << "num_unknowns=" << num_unknowns
                         << " max_mpi_message_size=" << max_message_bytes << ": " << blocks.size()
                         << " blocks, first of size " << blocks.front().size << ", last of size "
                         << blocks.back().size;
  EXPECT_LE(max_bytes, max_message_bytes)
    << "num_unknowns=" << num_unknowns << " max_mpi_message_size=" << max_message_bytes;
  EXPECT_TRUE(contiguous) << "num_unknowns=" << num_unknowns
                          << " max_mpi_message_size=" << max_message_bytes;
  EXPECT_TRUE(inside_buffer) << "num_unknowns=" << num_unknowns
                             << " max_mpi_message_size=" << max_message_bytes;
  EXPECT_EQ(next_position, static_cast<long long>(num_unknowns))
    << "max_mpi_message_size=" << max_message_bytes;
  EXPECT_EQ(msg_received.front().size(), blocks.size());
  EXPECT_EQ(max_num_messages, static_cast<int>(blocks.size()));
}

/// Calls SetupMessageData for one location with `num_unknowns` doubles, as for a successor
/// location (no receive flags), and returns its blocks. `max_num_messages` receives the count.
std::vector<AAH_MessageDetails>
SplitOneLocation(const SweepCommunicator& sweep_communicator,
                 const std::size_t num_unknowns,
                 const int max_message_bytes,
                 int& max_num_messages)
{
  const std::vector<int> locations = {mpi_comm.rank()};
  std::vector<std::vector<AAH_MessageDetails>> msg_data;
  max_num_messages = 0;
  SetupMessageData(
    locations,
    [num_unknowns](std::size_t) { return num_unknowns; },
    msg_data,
    nullptr,
    max_num_messages,
    sweep_communicator,
    max_message_bytes);
  if (msg_data.size() != 1)
    return {};
  return msg_data.front();
}

/// A message block: its size and its position in the buffer, in doubles.
struct ExpectedBlock
{
  long long size;
  long long position;
};

/// The split taken from the byte count: ceil(8 n / B) messages (one if the data fit in one
/// message, at most n), each of ceil(n / count) doubles but the last, which takes the rest. For a
/// limit B that is a whole number of doubles these blocks are valid.
std::vector<ExpectedBlock>
ByteBasedSplit(const long long num_unknowns, const long long max_message_bytes)
{
  const long long total_bytes = num_unknowns * static_cast<long long>(sizeof(double));
  long long count = 1;
  if (total_bytes > max_message_bytes)
    count = (total_bytes + max_message_bytes - 1) / max_message_bytes;
  count = std::min(count, num_unknowns);
  const long long size = (num_unknowns + count - 1) / count;

  std::vector<ExpectedBlock> blocks;
  blocks.reserve(static_cast<std::size_t>(count));
  for (long long m = 0; m + 1 < count; ++m)
    blocks.push_back({size, m * size});
  blocks.push_back({num_unknowns - ((count - 1) * size), (count - 1) * size});
  return blocks;
}

/// Equal blocks of `size` doubles, back to back from position 0.
std::vector<ExpectedBlock>
EqualBlocks(const long long count, const long long size)
{
  std::vector<ExpectedBlock> blocks;
  blocks.reserve(static_cast<std::size_t>(count));
  for (long long m = 0; m < count; ++m)
    blocks.push_back({size, m * size});
  return blocks;
}

/// Expects `blocks` to equal `expected`, and reports the first block that differs.
void
ExpectBlocks(const std::vector<AAH_MessageDetails>& blocks,
             const std::vector<ExpectedBlock>& expected,
             const std::size_t num_unknowns,
             const int max_message_bytes)
{
  ASSERT_EQ(blocks.size(), expected.size())
    << "num_unknowns=" << num_unknowns << " max_mpi_message_size=" << max_message_bytes;
  for (std::size_t m = 0; m < blocks.size(); ++m)
  {
    if (blocks[m].size != expected[m].size or blocks[m].block_pos != expected[m].position)
    {
      ADD_FAILURE() << "num_unknowns=" << num_unknowns
                    << " max_mpi_message_size=" << max_message_bytes << ": block " << m << " has "
                    << blocks[m].size << " doubles at " << blocks[m].block_pos << ", expected "
                    << expected[m].size << " at " << expected[m].position;
      return;
    }
  }
}

} // namespace

TEST(AAHMessageStructTest, EmptyLocationHasNoMessages)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);
  const std::vector<int> locations = {mpi_comm.rank()};
  std::vector<std::vector<AAH_MessageDetails>> msg_data;
  std::vector<std::vector<bool>> msg_received;
  int max_num_messages = 0;
  SetupMessageData(
    locations,
    [](std::size_t) { return std::size_t{0}; },
    msg_data,
    &msg_received,
    max_num_messages,
    sweep_communicator,
    32768);

  ASSERT_EQ(msg_data.size(), 1U);
  EXPECT_TRUE(msg_data.front().empty());
  EXPECT_TRUE(msg_received.front().empty());
  EXPECT_EQ(max_num_messages, 0);
}

TEST(AAHMessageStructTest, MessageSizeNotMultipleOfDouble)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);

  // max_mpi_message_size accepts any value >= 1024 bytes. With 1025 bytes the message count
  // ceil(8 n / 1025) and the block size ceil(n / count) = 129 doubles (1032 bytes) leave a last
  // block of 0 doubles for n = 16641, of -1 for n = 16769 and of -620 for n = 100000, whose
  // preceding blocks also extend past the end of the buffer.
  ExpectValidSplit(sweep_communicator, 16641, 1025);
  ExpectValidSplit(sweep_communicator, 16769, 1025);
  ExpectValidSplit(sweep_communicator, 100000, 1025);
}

TEST(AAHMessageStructTest, BlocksPartitionUnknowns)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);

  for (const int max_message_bytes : {1024, 1025, 1028, 1031, 1032, 32768, 32769})
    for (const std::size_t num_unknowns : {1, 2, 127, 128, 129, 1000, 4096, 16641, 16769, 100000})
      ExpectValidSplit(sweep_communicator, num_unknowns, max_message_bytes);
}

TEST(AAHMessageStructTest, WholeDoubleLimitsKeepByteBasedSplit)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);

  // A limit that is a whole number of doubles, the default 32768 bytes included, keeps the split
  // taken from the byte count: the same messages, block sizes and positions.
  for (const int max_message_bytes : {1024, 1032, 4096, 32768, 262144, 2147483640})
    for (const std::size_t num_unknowns :
         {1, 2, 127, 128, 129, 1000, 4095, 4096, 4097, 8193, 16641, 16769, 18704, 100000, 1000000})
    {
      int max_num_messages = 0;
      const auto blocks =
        SplitOneLocation(sweep_communicator, num_unknowns, max_message_bytes, max_num_messages);
      const auto expected = ByteBasedSplit(static_cast<long long>(num_unknowns), max_message_bytes);
      ExpectBlocks(blocks, expected, num_unknowns, max_message_bytes);
      EXPECT_EQ(max_num_messages, static_cast<int>(expected.size()))
        << "num_unknowns=" << num_unknowns << " max_mpi_message_size=" << max_message_bytes;
    }

  // At 32768 bytes (4096 doubles), 18704 doubles go in five messages of 3741, 3741, 3741, 3741
  // and 3740 doubles, and 100000 doubles in 25 messages of 4000 doubles.
  int max_num_messages = 0;
  ExpectBlocks(SplitOneLocation(sweep_communicator, 18704, 32768, max_num_messages),
               {{3741, 0}, {3741, 3741}, {3741, 7482}, {3741, 11223}, {3740, 14964}},
               18704,
               32768);
  ExpectBlocks(SplitOneLocation(sweep_communicator, 100000, 32768, max_num_messages),
               EqualBlocks(25, 4000),
               100000,
               32768);
}

TEST(AAHMessageStructTest, MessageCountIsMinimal)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);

  // A block holds at most floor(B / 8) whole doubles, so no valid split of n doubles has fewer
  // than ceil(n / floor(B / 8)) messages. The split uses that many, whether or not B is a whole
  // number of doubles.
  for (const int max_message_bytes :
       {1024, 1025, 1028, 1031, 1032, 1033, 4100, 32767, 32768, 32769, 32775, 2147483647})
    for (const std::size_t num_unknowns :
         {1, 2, 127, 128, 129, 1000, 4095, 4096, 4097, 8193, 16641, 16769, 18704, 100000, 1000000})
    {
      ExpectValidSplit(sweep_communicator, num_unknowns, max_message_bytes);
      int max_num_messages = 0;
      const auto blocks =
        SplitOneLocation(sweep_communicator, num_unknowns, max_message_bytes, max_num_messages);
      const std::size_t max_doubles = static_cast<std::size_t>(max_message_bytes) / sizeof(double);
      EXPECT_EQ(blocks.size(), (num_unknowns + max_doubles - 1) / max_doubles)
        << "num_unknowns=" << num_unknowns << " max_mpi_message_size=" << max_message_bytes;
    }
}

TEST(AAHMessageStructTest, LargestLimitsWithLargeLocation)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);

  // At the largest int limits, 2^31 - 8 and 2^31 - 1 bytes (268435455 doubles), two billion
  // doubles go in eight messages of 250000000 doubles; the count must not overflow int.
  int max_num_messages = 0;
  ExpectBlocks(SplitOneLocation(sweep_communicator, 2000000000, 2147483640, max_num_messages),
               EqualBlocks(8, 250000000),
               2000000000,
               2147483640);
  EXPECT_EQ(max_num_messages, 8);

  ExpectValidSplit(sweep_communicator, 2000000000, 2147483647);
  EXPECT_EQ(SplitOneLocation(sweep_communicator, 2000000000, 2147483647, max_num_messages).size(),
            8U);
}

TEST(AAHMessageStructTest, UnknownCountsUpToIntMax)
{
  SweepCommunicator sweep_communicator(mpi_comm, 1);

  // A location may hold up to INT_MAX = 2147483647 unknowns; one more is refused. Near that count
  // the rounded-up sum n + count - 1 of a ceiling division exceeds INT_MAX. Each case gives n, the
  // limit B and the minimal message count ceil(n / floor(B / 8)). That sum first exceeds INT_MAX
  // at n = 2147483641 for B = 2^31 - 8 and 2^31 - 1 (268435455 doubles per message), at
  // n = 2147483633 for B = 2^30 and 2^30 + 7 (2^27 doubles) and at n = 2147482625 for B = 2^24
  // and 2^24 + 7 (2^21 doubles). The count one below each of these is included; 2147483640
  // doubles fill eight messages of 2^31 - 8 bytes exactly.
  struct Case
  {
    std::size_t num_unknowns;
    int max_message_bytes;
    std::size_t num_messages;
  };
  const std::vector<Case> cases = {{2147483647, 2147483647, 9},
                                   {2147483647, 2147483640, 9},
                                   {2147483641, 2147483647, 9},
                                   {2147483640, 2147483640, 8},
                                   {2147483647, 1073741824, 16},
                                   {2147483633, 1073741831, 16},
                                   {2147483632, 1073741824, 16},
                                   {2147483647, 16777216, 1024},
                                   {2147482625, 16777223, 1024},
                                   {2147482624, 16777216, 1024}};
  int max_num_messages = 0;
  for (const auto& c : cases)
  {
    ExpectValidSplit(sweep_communicator, c.num_unknowns, c.max_message_bytes);
    EXPECT_EQ(
      SplitOneLocation(sweep_communicator, c.num_unknowns, c.max_message_bytes, max_num_messages)
        .size(),
      c.num_messages)
      << "num_unknowns=" << c.num_unknowns << " max_mpi_message_size=" << c.max_message_bytes;
  }

  EXPECT_THROW(SplitOneLocation(sweep_communicator, 2147483648, 2147483647, max_num_messages),
               std::logic_error);
}
