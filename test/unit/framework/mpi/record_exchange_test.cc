// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mpi/record_exchange.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

/// Payload of the record with `key` sent by `source`: `key % 4` values derived from both.
std::vector<double>
MakePayload(const std::uint64_t key, const int source)
{
  std::vector<double> payload;
  payload.reserve(key % 4);
  for (std::uint64_t i = 0; i < key % 4; ++i)
    payload.push_back(static_cast<double>(1000 * source) + static_cast<double>(key) +
                      0.5 * static_cast<double>(i));
  return payload;
}

/// Where a FakeFileSet fails, on its failing rank.
enum class FailAt
{
  never,
  open,
  prepare,
  begin,
  fill,
  receive,
  /// FillChunk addresses a record to a rank that does not exist, so ExchangeRecords rejects it.
  owner,
  /// PrepareFile throws std::logic_error on every rank, as a rejected owner lookup does.
  prepare_throw,
  header_width,
  done
};

/// What a FakeFileSet throws where it fails.
enum class FailWith
{
  standard,
  non_standard
};

/// A thrown object that is not a std::exception.
struct NonStandardFailure
{
};

/**
 * An in-memory file set for RunFileSetPasses: file k holds records with keys 10 k + i for
 * i < k % 3 (so some files are empty), sent in two streams to rank Owner(key, P). Stream s
 * carries the payload {key + s}.
 */
class FakeFileSet
{
public:
  FakeFileSet(const FailAt fail_at,
              const int failing_rank,
              const FailWith fail_with = FailWith::standard)
    : fail_at_(fail_at), failing_rank_(failing_rank), fail_with_(fail_with)
  {
  }

  /// Record 10 k + i goes to rank (k + i + 1) mod P: with 3 P files, every rank, including the
  /// failing rank P - 1, receives records for every P from 1 to 64. For P = 1 and 2 some records
  /// stay on the rank that reads their file, which also exercises sends to self.
  static std::uint64_t Owner(const std::uint64_t key, const std::uint64_t num_ranks)
  {
    return (key / 10 + key % 10 + 1) % num_ranks;
  }

  static std::vector<std::uint64_t> FileKeys(const std::uint64_t file_id)
  {
    std::vector<std::uint64_t> keys;
    keys.reserve(file_id % 3);
    for (std::uint64_t i = 0; i < file_id % 3; ++i)
      keys.push_back(10 * file_id + i);
    return keys;
  }

  void OpenFile(const std::uint64_t file_id)
  {
    ++num_opened_;
    Fail(FailAt::open);
    keys_ = FileKeys(file_id);
  }

  std::string PrepareFile(const bool has_file)
  {
    // A collective step on every rank, as a real source's owner lookup is.
    int files = has_file ? 1 : 0;
    mpi_comm.all_reduce(files, mpi::op::sum<int>());
    if (fail_at_ == FailAt::prepare_throw)
      throw std::logic_error("injected collective prepare failure");
    if (fail_at_ == FailAt::prepare and mpi_comm.rank() == failing_rank_)
      return "injected prepare failure";
    return {};
  }

  std::size_t GetHeaderWidth(std::size_t /*stream*/) const
  {
    Fail(FailAt::header_width);
    return 1;
  }

  void BeginStream(std::size_t /*stream*/)
  {
    Fail(FailAt::begin);
    next_ = 0;
  }

  bool IsStreamDone(std::size_t /*stream*/) const
  {
    Fail(FailAt::done);
    return next_ == keys_.size();
  }

  void FillChunk(const std::size_t stream, std::uint64_t /*budget*/, RecordBatch& batch)
  {
    Fail(FailAt::fill);
    // One record per chunk, the smallest the budget rule allows.
    const auto key = keys_[next_++];
    const auto num_ranks = static_cast<std::uint64_t>(mpi_comm.size());
    const bool bad_owner = fail_at_ == FailAt::owner and mpi_comm.rank() == failing_rank_;
    batch.AddRecord(static_cast<int>(bad_owner ? num_ranks : Owner(key, num_ranks)))[0] =
      static_cast<std::int64_t>(key);
    const auto value = static_cast<double>(key + stream);
    batch.AppendPayload(&value, 1);
  }

  void Receive(const std::size_t stream, const std::uint64_t pass, const ReceivedRecords& records)
  {
    if (records.GetNumRecords() > 0)
      Fail(FailAt::receive);
    for (std::size_t r = 0; r < records.GetNumRecords(); ++r)
    {
      const auto key = static_cast<std::uint64_t>(records.GetHeader(r)[0]);
      const auto file_id = pass * static_cast<std::uint64_t>(mpi_comm.size()) +
                           static_cast<std::uint64_t>(records.sources[r]);
      received_.emplace_back(stream, key, file_id, records.GetPayload(r)[0]);
    }
  }

  void CloseFile()
  {
    ++num_closed_;
    keys_.clear();
  }

  /// OpenFile and CloseFile calls on this rank.
  std::uint64_t GetNumOpened() const { return num_opened_; }
  std::uint64_t GetNumClosed() const { return num_closed_; }

  /// (stream, key, file, payload) of every record received.
  const std::vector<std::tuple<std::size_t, std::uint64_t, std::uint64_t, double>>&
  GetReceived() const
  {
    return received_;
  }

private:
  void Fail(const FailAt where) const
  {
    if (fail_at_ != where or mpi_comm.rank() != failing_rank_)
      return;
    if (fail_with_ == FailWith::non_standard)
      throw NonStandardFailure{};
    throw std::runtime_error("injected failure");
  }

  const FailAt fail_at_;
  const int failing_rank_;
  const FailWith fail_with_;
  std::vector<std::uint64_t> keys_;
  std::size_t next_ = 0;
  std::uint64_t num_opened_ = 0;
  std::uint64_t num_closed_ = 0;
  std::vector<std::tuple<std::size_t, std::uint64_t, std::uint64_t, double>> received_;
};

} // namespace

TEST(RecordExchangeTest, DeliversRecordsToTheirOwners)
{
  const auto rank = mpi_comm.rank();
  const auto num_ranks = mpi_comm.size();

  // Rank r sends keys r, r + P, ..., some to itself, with payloads of length 0 to 3. With more than
  // one rank, rank 0 sends nothing, so a rank may have no records; alone, it sends to itself.
  RecordBatch batch(2);
  std::vector<std::uint64_t> sent_keys;
  if (rank > 0 or num_ranks == 1)
    for (auto key = static_cast<std::uint64_t>(rank); key < 40;
         key += static_cast<std::uint64_t>(num_ranks))
      sent_keys.push_back(key);
  for (const auto key : sent_keys)
  {
    auto* header = batch.AddRecord(static_cast<int>((key * 7) % num_ranks));
    header[0] = static_cast<std::int64_t>(key);
    header[1] = -static_cast<std::int64_t>(key);
    const auto payload = MakePayload(key, rank);
    batch.AppendPayload(payload.data(), payload.size());
  }
  const auto received = ExchangeRecords(mpi_comm, std::move(batch), "Test exchange");

  // Every key with (7 key) mod P == rank arrives once, from the rank that sent it.
  std::set<std::uint64_t> expected;
  for (int source = num_ranks == 1 ? 0 : 1; source < num_ranks; ++source)
    for (auto key = static_cast<std::uint64_t>(source); key < 40;
         key += static_cast<std::uint64_t>(num_ranks))
      if (std::cmp_equal((key * 7) % num_ranks, rank))
        expected.insert(key);
  ASSERT_EQ(received.GetNumRecords(), expected.size());
  std::set<std::uint64_t> seen;
  int previous_source = 0;
  for (std::size_t r = 0; r < received.GetNumRecords(); ++r)
  {
    const auto* header = received.GetHeader(r);
    const auto key = static_cast<std::uint64_t>(header[0]);
    EXPECT_EQ(header[1], -header[0]);
    EXPECT_EQ(received.sources[r], static_cast<int>(key % num_ranks));
    // Arrival order is by source rank.
    EXPECT_GE(received.sources[r], previous_source);
    previous_source = received.sources[r];
    const auto payload = MakePayload(key, received.sources[r]);
    ASSERT_EQ(received.GetPayloadSize(r), payload.size());
    for (std::size_t i = 0; i < payload.size(); ++i)
      EXPECT_EQ(received.GetPayload(r)[i], payload[i]);
    EXPECT_TRUE(seen.insert(key).second);
  }
  EXPECT_EQ(seen, expected);
}

// Each test below makes one exchange and checks nothing until after it, so a failed fatal
// assertion cannot leave the other ranks waiting in a later collective.
TEST(RecordExchangeTest, DeliversRecordsWithoutHeaderValues)
{
  const auto rank = mpi_comm.rank();
  const auto num_ranks = mpi_comm.size();

  // Records without header values, each with the payload {rank, destination}.
  RecordBatch payload_only(0);
  for (int destination = 0; destination < num_ranks; ++destination)
  {
    payload_only.AddRecord(destination);
    const std::array<double, 2> payload = {static_cast<double>(rank),
                                           static_cast<double>(destination)};
    payload_only.AppendPayload(payload.data(), payload.size());
  }
  const auto received = ExchangeRecords(mpi_comm, std::move(payload_only), "Test exchange");
  ASSERT_EQ(received.GetNumRecords(), static_cast<std::size_t>(num_ranks));
  // One record from each rank, carrying that rank's payload.
  std::multiset<int> sources(received.sources.begin(), received.sources.end());
  std::multiset<int> all_ranks;
  for (int p = 0; p < num_ranks; ++p)
    all_ranks.insert(p);
  EXPECT_EQ(sources, all_ranks);
  for (std::size_t r = 0; r < received.GetNumRecords(); ++r)
  {
    ASSERT_EQ(received.GetPayloadSize(r), 2U);
    EXPECT_EQ(received.GetPayload(r)[0], static_cast<double>(received.sources[r]));
    EXPECT_EQ(received.GetPayload(r)[1], static_cast<double>(rank));
  }
}

TEST(RecordExchangeTest, DeliversRecordsWithoutPayload)
{
  const auto rank = mpi_comm.rank();
  const auto num_ranks = mpi_comm.size();

  // Records without payload, so no rank sends or receives any payload value.
  RecordBatch header_only(1);
  for (int destination = 0; destination < num_ranks; ++destination)
    header_only.AddRecord(destination)[0] = rank;
  const auto headers = ExchangeRecords(mpi_comm, std::move(header_only), "Test exchange");
  ASSERT_EQ(headers.GetNumRecords(), static_cast<std::size_t>(num_ranks));
  EXPECT_TRUE(headers.payload.empty());
  // One record from each rank, carrying that rank's id.
  std::multiset<std::int64_t> ids;
  std::multiset<std::int64_t> all_ranks;
  for (std::size_t r = 0; r < headers.GetNumRecords(); ++r)
  {
    EXPECT_EQ(headers.GetHeader(r)[0], headers.sources[r]);
    EXPECT_EQ(headers.GetPayloadSize(r), 0U);
    ids.insert(headers.GetHeader(r)[0]);
  }
  for (int p = 0; p < num_ranks; ++p)
    all_ranks.insert(p);
  EXPECT_EQ(ids, all_ranks);
}

TEST(RecordExchangeTest, RejectsAHeaderWidthWithoutRoomForThePayloadLength)
{
  // The packed width adds one value to the header width, which must not wrap to zero, even for an
  // empty batch.
  RecordBatch batch(std::numeric_limits<std::size_t>::max());
  EXPECT_THROW(ExchangeRecords(mpi_comm, std::move(batch), "Test exchange"), std::logic_error);
}

TEST(RecordExchangeTest, RejectsAnOwnerOutsideTheCommunicator)
{
  RecordBatch batch(1);
  if (mpi_comm.rank() == 0)
    batch.AddRecord(mpi_comm.size())[0] = 1;
  EXPECT_THROW(ExchangeRecords(mpi_comm, std::move(batch), "Test exchange"), std::logic_error);
}

TEST(RecordExchangeTest, OwnerDirectoryFindsRegisteredKeys)
{
  const auto rank = mpi_comm.rank();
  const auto num_ranks = static_cast<std::uint64_t>(mpi_comm.size());

  // Rank r owns keys 3 k + 1 with k mod P == r; every rank asks for keys 0 to 3 * 10 P.
  std::vector<std::uint64_t> owned;
  for (auto k = static_cast<std::uint64_t>(rank); k < 10 * num_ranks; k += num_ranks)
    owned.push_back(3 * k + 1);
  const OwnerDirectory directory(mpi_comm, owned);

  std::vector<std::uint64_t> keys;
  keys.reserve(30 * num_ranks);
  for (std::uint64_t key = 0; key < 30 * num_ranks; ++key)
    keys.push_back(key);
  const auto owners = directory.FindOwners(keys);
  // A nonfatal check, so that every rank reaches the second lookup below.
  EXPECT_EQ(owners.size(), keys.size());
  for (std::size_t i = 0; i < std::min(keys.size(), owners.size()); ++i)
  {
    const auto expected = keys[i] % 3 == 1 ? static_cast<int>(((keys[i] - 1) / 3) % num_ranks) : -1;
    EXPECT_EQ(owners[i], expected) << "key " << keys[i];
  }

  // A rank may ask for nothing while the others ask.
  EXPECT_TRUE(directory.FindOwners(rank == 0 ? std::vector<std::uint64_t>{} : keys).size() ==
              (rank == 0 ? 0 : keys.size()));
}

TEST(RecordExchangeTest, OwnerDirectoryRejectsDuplicateKeys)
{
  // Every rank registers key 5, and rank 0 registers it twice in a one-rank run.
  std::vector<std::uint64_t> owned = {5};
  if (mpi_comm.size() == 1)
    owned.push_back(5);
  EXPECT_THROW(OwnerDirectory(mpi_comm, owned), std::logic_error);
}

TEST(RecordExchangeTest, FileSetPassesDeliverEveryRecordOnce)
{
  const auto rank = static_cast<std::uint64_t>(mpi_comm.rank());
  const auto num_ranks = static_cast<std::uint64_t>(mpi_comm.size());
  for (const auto num_files : {std::uint64_t{1}, num_ranks + 1, 2 * num_ranks + 1, num_ranks})
  {
    FakeFileSet source(FailAt::never, 0);
    RunFileSetPasses<std::runtime_error>(mpi_comm, num_files, 2, 1, source, "Test file set");

    std::multiset<std::tuple<std::size_t, std::uint64_t, std::uint64_t, double>> expected, seen;
    for (std::uint64_t file_id = 0; file_id < num_files; ++file_id)
      for (const auto key : FakeFileSet::FileKeys(file_id))
        if (FakeFileSet::Owner(key, num_ranks) == rank)
          for (std::size_t stream = 0; stream < 2; ++stream)
            expected.emplace(stream, key, file_id, static_cast<double>(key + stream));
    for (const auto& record : source.GetReceived())
      seen.insert(record);
    EXPECT_EQ(seen, expected) << num_files << " files";
    // Each rank opens exactly the files k = rank (mod P), including files without records.
    std::uint64_t expected_opens = 0;
    for (std::uint64_t file_id = 0; file_id < num_files; ++file_id)
      if (file_id % num_ranks == rank)
        ++expected_opens;
    EXPECT_EQ(source.GetNumOpened(), expected_opens) << num_files << " files";
    EXPECT_EQ(source.GetNumClosed(), source.GetNumOpened()) << num_files << " files";
  }
}

TEST(RecordExchangeTest, FileSetPassesRaiseAFailureOnEveryRank)
{
  // Rank P - 1 fails at each step in turn. With 3 P files, rank P - 1 always reads a file that has
  // records (file 3 P - 1), and every rank receives a record (see FakeFileSet::Owner).
  const auto num_files = 3 * static_cast<std::uint64_t>(mpi_comm.size());
  // FailAt::owner makes ExchangeRecords throw std::logic_error on every rank, which is raised as
  // the requested exception type.
  // FailAt::prepare_throw does the same from PrepareFile.
  for (const auto fail_at : {FailAt::open,
                             FailAt::prepare,
                             FailAt::prepare_throw,
                             FailAt::header_width,
                             FailAt::begin,
                             FailAt::done,
                             FailAt::fill,
                             FailAt::receive,
                             FailAt::owner})
  {
    FakeFileSet source(fail_at, mpi_comm.size() - 1);
    EXPECT_THROW(
      RunFileSetPasses<std::runtime_error>(mpi_comm, num_files, 2, 1, source, "Test file set"),
      std::runtime_error);
    // A file that was opened is closed, whichever step failed.
    EXPECT_EQ(source.GetNumClosed(), source.GetNumOpened());
    // Every rank left the passes, so a following collective completes.
    int count = 1;
    mpi_comm.all_reduce(count, mpi::op::sum<int>());
    EXPECT_EQ(count, mpi_comm.size());
  }
}

TEST(RecordExchangeTest, FileSetPassesRaiseANonStandardFailureOnEveryRank)
{
  // As above, but rank P - 1 throws an object that is not a std::exception.
  const auto num_files = 3 * static_cast<std::uint64_t>(mpi_comm.size());
  for (const auto fail_at : {FailAt::open,
                             FailAt::header_width,
                             FailAt::begin,
                             FailAt::done,
                             FailAt::fill,
                             FailAt::receive})
  {
    FakeFileSet source(fail_at, mpi_comm.size() - 1, FailWith::non_standard);
    EXPECT_THROW(
      RunFileSetPasses<std::runtime_error>(mpi_comm, num_files, 2, 1, source, "Test file set"),
      std::runtime_error);
    EXPECT_EQ(source.GetNumClosed(), source.GetNumOpened());
    int count = 1;
    mpi_comm.all_reduce(count, mpi::op::sum<int>());
    EXPECT_EQ(count, mpi_comm.size());
  }
}
