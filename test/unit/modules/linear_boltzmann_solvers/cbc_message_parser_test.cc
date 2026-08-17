// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_message_parser.h"
#include "gtest/gtest.h"
#include <cstring>
#include <limits>
#include <stdexcept>

using namespace opensn;

namespace
{
template <typename T>
void
Put(std::vector<char>& packet, T value)
{
  const auto size = packet.size();
  packet.resize(size + sizeof(T));
  std::memcpy(packet.data() + size, &value, sizeof(T));
}

std::vector<char>
Record(std::uint8_t kind,
       std::size_t slot,
       std::size_t total = 3,
       std::size_t offset = 0,
       std::size_t count = 3)
{
  std::vector<char> packet;
  Put(packet, kind);
  Put(packet, slot);
  if (kind == 2 or kind == 3)
  {
    Put(packet, total);
    Put(packet, offset);
    Put(packet, count);
  }
  for (std::size_t i = 0; i < count; ++i)
    Put(packet, static_cast<double>(i + offset + 1));
  return packet;
}

detail::CBCMessageParser
Parser()
{
  return detail::CBCMessageParser({{0, 3}}, {{1, 3}}, std::vector<int>{1}, 2);
}
} // namespace

TEST(CBCMessageParser, WholeFragmentsCompletionAndReset)
{
  auto parser = Parser();
  std::vector<double> values(3);
  int completions = 0;
  const auto Store = [&](const detail::CBCMessageParser::Record& r)
  {
    std::memcpy(values.data() + r.offset, r.payload, r.count * sizeof(double));
    completions += r.complete;
  };
  for (int repeat = 0; repeat < 2; ++repeat)
  {
    parser.Reset();
    EXPECT_NO_THROW(parser.Parse(Record(0, 0), 0, false, Store));
    EXPECT_NO_THROW(parser.Parse(Record(3, 0, 3, 0, 1), 1, true, Store));
    EXPECT_NO_THROW(parser.Parse(Record(3, 0, 3, 1, 2), 1, true, Store));
    EXPECT_EQ(values, (std::vector<double>{1, 2, 3}));
    EXPECT_FALSE(parser.DelayedComplete());
    EXPECT_NO_THROW(parser.Parse(std::vector<char>{4}, 1, true, Store));
    EXPECT_TRUE(parser.DelayedComplete());
  }
  EXPECT_EQ(completions, 4);
}

TEST(CBCMessageParser, EveryTruncatedHeaderAndPayload)
{
  const auto Store = [](const auto&) {};
  for (const std::uint8_t kind : {0, 1, 2, 3})
  {
    auto packet = Record(kind, 0);
    for (std::size_t length = 0; length < packet.size(); ++length)
    {
      auto parser = Parser();
      EXPECT_THROW(parser.Parse(std::span(packet).first(length), kind % 2, kind % 2, Store),
                   std::logic_error)
        << static_cast<int>(kind) << " length " << length;
    }
  }
}

TEST(CBCMessageParser, InvalidMetadataKindsSourcesSlotsAndBounds)
{
  const auto Store = [](const auto&) {};
  EXPECT_THROW((detail::CBCMessageParser({{-1, 3}}, {}, {}, 2)), std::logic_error);
  EXPECT_THROW((detail::CBCMessageParser({}, {{1, 3}}, {}, 2)), std::logic_error);
  for (const auto& packet :
       std::vector<std::vector<char>>{{5},
                                      Record(0, 1),
                                      Record(2, 0, 2),
                                      Record(2, 0, 3, 1, 2),
                                      Record(2, 0, 3, 0, 0),
                                      Record(2, 0, std::numeric_limits<std::size_t>::max()),
                                      Record(2, 0, 3, 0, 4)})
  {
    auto parser = Parser();
    EXPECT_THROW(parser.Parse(packet, 0, false, Store), std::logic_error);
  }
  auto parser = Parser();
  EXPECT_THROW(parser.Parse(Record(0, 0), -1, false, Store), std::logic_error);
  EXPECT_THROW(parser.Parse(Record(0, 0), 2, false, Store), std::logic_error);
  EXPECT_THROW(parser.Parse(Record(0, 0), 1, false, Store), std::logic_error);
  EXPECT_THROW(parser.Parse(Record(1, 0), 1, false, Store), std::logic_error);
}

TEST(CBCMessageParser, DuplicateOverlapAndCompletionMisuse)
{
  const auto Store = [](const auto&) {};
  auto parser = Parser();
  EXPECT_THROW(parser.Parse(std::vector<char>{4}, 1, true, Store), std::logic_error);
  EXPECT_THROW(parser.Parse(std::vector<char>{4}, 0, false, Store), std::logic_error);
  EXPECT_NO_THROW(parser.Parse(Record(2, 0, 3, 0, 2), 0, false, Store));
  EXPECT_THROW(parser.Parse(Record(2, 0, 3, 0, 2), 0, false, Store), std::logic_error);
  EXPECT_THROW(parser.Parse(Record(2, 0, 3, 1, 2), 0, false, Store), std::logic_error);
  EXPECT_NO_THROW(parser.Parse(Record(2, 0, 3, 2, 1), 0, false, Store));
  EXPECT_THROW(parser.Parse(Record(0, 0), 0, false, Store), std::logic_error);
  EXPECT_NO_THROW(parser.Parse(Record(1, 0), 1, true, Store));
  EXPECT_NO_THROW(parser.Parse(std::vector<char>{4}, 1, true, Store));
  EXPECT_THROW(parser.Parse(std::vector<char>{4}, 1, true, Store), std::logic_error);
  EXPECT_THROW(parser.Parse(Record(1, 0), 1, true, Store), std::logic_error);
}
