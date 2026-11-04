// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/runtime.h"
#include "framework/utils/collective_error.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace opensn
{

/**
 * Records to send to their owner ranks with ExchangeRecords. Each record has an owner rank,
 * `header_width` int64 header values, and a payload of any number of doubles.
 */
struct RecordBatch
{
  explicit RecordBatch(const std::size_t header_width) : header_width(header_width) {}

  /// Append a record for rank `owner` with an empty payload. Returns its header to fill.
  std::int64_t* AddRecord(int owner);
  /// Append `count` values to the payload of the last record.
  void AppendPayload(const double* values, std::size_t count);
  std::size_t GetNumRecords() const { return owners.size(); }
  /// Bytes that the records occupy in an exchange.
  std::uint64_t GetPackedBytes() const;
  /// Bytes that one record with `payload_size` payload values occupies in an exchange.
  std::uint64_t GetRecordBytes(std::uint64_t payload_size) const;

  std::size_t header_width;
  std::vector<int> owners;
  std::vector<std::int64_t> headers;
  std::vector<std::uint64_t> payload_offsets = {0};
  std::vector<double> payload;
};

/// Records received by ExchangeRecords, ordered by source rank, then in the order they were sent.
struct ReceivedRecords
{
  std::size_t header_width = 0;
  /// Source rank of each record.
  std::vector<int> sources;
  std::vector<std::int64_t> headers;
  std::vector<std::uint64_t> payload_offsets = {0};
  std::vector<double> payload;

  std::size_t GetNumRecords() const { return sources.size(); }
  const std::int64_t* GetHeader(const std::size_t record) const
  {
    return headers.data() + record * header_width;
  }
  const double* GetPayload(const std::size_t record) const
  {
    return payload.data() + payload_offsets[record];
  }
  std::uint64_t GetPayloadSize(const std::size_t record) const
  {
    return payload_offsets[record + 1] - payload_offsets[record];
  }
};

/**
 * Send every record of `records` to its owner rank and return the records this rank receives.
 * Collective over `comm`; every rank's batch must have the same header width. Each record's
 * payload length travels with it, so payloads may differ in length. Throws std::logic_error on
 * every rank, before any data are sent, if an owner is not a rank of `comm` or if the values a rank
 * sends or receives exceed the MPI int count range, and after the exchange if the received data are
 * inconsistent. `operation` describes the caller in error messages.
 */
ReceivedRecords
ExchangeRecords(const mpi::Communicator& comm, RecordBatch&& records, const std::string& operation);

/**
 * A distributed map from keys to the rank that owns each key. Key k is kept by rank k mod P, so
 * each rank stores the entries of about 1/P of the keys and no rank holds a global table.
 */
class OwnerDirectory
{
public:
  /**
   * Register `local_keys`, the keys owned by this rank. Collective over `comm`. Throws
   * std::logic_error on every rank if a key is registered more than once.
   */
  OwnerDirectory(const mpi::Communicator& comm, const std::vector<std::uint64_t>& local_keys);

  /// Return the owner of each of `keys`, or -1 for a key that no rank registered. Collective.
  std::vector<int> FindOwners(const std::vector<std::uint64_t>& keys) const;

private:
  mpi::Communicator comm_;
  /// (key, owner) pairs of the keys this rank keeps, sorted by key.
  std::vector<std::pair<std::uint64_t, int>> entries_;
};

/**
 * Read a set of `num_files` files, one per writing rank, on the ranks of `comm` and send their
 * records to their owners. File k is read by rank k mod P, in pass k / P, so each rank reads at
 * most ceil(num_files / P) files, one at a time. Each file is sent as `num_streams` streams of
 * records in chunks of at most `max_send_bytes` packed bytes (at least one record). Every rank runs
 * every pass and every exchange, including ranks without a file.
 *
 * `source` provides:
 * - `void OpenFile(std::uint64_t file_id)`: open file `file_id`. Rank-local.
 * - `std::string PrepareFile(bool has_file)`: called on every rank after the open. Collective;
 *   returns this rank's error, if any.
 * - `std::size_t GetHeaderWidth(std::size_t stream) const`: the same on every rank. Must not
 *   throw.
 * - `void BeginStream(std::size_t stream)`: called with an open file. Rank-local.
 * - `bool IsStreamDone(std::size_t stream) const`: called with an open file. Must not throw.
 * - `void FillChunk(std::size_t stream, std::uint64_t budget, RecordBatch& batch)`: add the next
 *   records that fit in `budget` bytes, at least one. Rank-local.
 * - `void Receive(std::size_t stream, std::uint64_t pass, const ReceivedRecords& records)`:
 *   take the records this rank owns; a record from rank r comes from file pass * P + r. Rank-local.
 * - `void CloseFile()`: release the open file. Called once after every OpenFile call, also when
 *   OpenFile or a later step throws, so it must accept a partly opened file. Rank-local; must not
 *   throw.
 *
 * The other rank-local callbacks may throw. Any error, including an exchange that ExchangeRecords
 * rejects, is raised as `Exception` on every rank, naming the lowest failing rank and `operation`,
 * so no rank is left waiting in a collective.
 */
template <typename Exception, typename Source>
void
RunFileSetPasses(const mpi::Communicator& comm,
                 const std::uint64_t num_files,
                 const std::size_t num_streams,
                 const std::uint64_t max_send_bytes,
                 Source& source,
                 const std::string& operation)
{
  const auto num_ranks = static_cast<std::uint64_t>(comm.size());
  const auto rank = static_cast<std::uint64_t>(comm.rank());
  const auto num_passes = (num_files + num_ranks - 1) / num_ranks;
  // Closes an opened file when the pass ends, whether it completes or throws.
  struct FileCloser
  {
    Source* source = nullptr;
    FileCloser() = default;
    FileCloser(const FileCloser&) = delete;
    FileCloser& operator=(const FileCloser&) = delete;
    ~FileCloser()
    {
      if (source != nullptr)
        source->CloseFile();
    }
  };

  for (std::uint64_t pass = 0; pass < num_passes; ++pass)
  {
    const auto file_id = pass * num_ranks + rank;
    const bool has_file = file_id < num_files;
    std::string error;
    FileCloser closer;
    if (has_file)
    {
      closer.source = &source;
      try
      {
        source.OpenFile(file_id);
      }
      catch (const std::exception& exception)
      {
        error = DescribeError(exception);
      }
    }
    ThrowIfAnyRankFailed<Exception>(comm, error, operation);
    error = source.PrepareFile(has_file);
    ThrowIfAnyRankFailed<Exception>(comm, error, operation);

    for (std::size_t stream = 0; stream < num_streams; ++stream)
    {
      if (has_file)
      {
        try
        {
          source.BeginStream(stream);
        }
        catch (const std::exception& exception)
        {
          error = DescribeError(exception);
        }
      }
      while (true)
      {
        RecordBatch batch(source.GetHeaderWidth(stream));
        if (error.empty() and has_file and not source.IsStreamDone(stream))
        {
          try
          {
            source.FillChunk(stream, max_send_bytes, batch);
          }
          catch (const std::exception& exception)
          {
            error = DescribeError(exception);
            batch = RecordBatch(batch.header_width);
          }
        }
        // Every rank takes part in every exchange until no rank has records left in the stream.
        std::array<int, 2> local_flags = {error.empty() ? 0 : 1,
                                          batch.GetNumRecords() == 0 ? 0 : 1};
        std::array<int, 2> global_flags = {0, 0};
        comm.all_reduce(local_flags.data(), 2, global_flags.data(), mpi::op::max<int>());
        if (global_flags[0] != 0)
          ThrowIfAnyRankFailed<Exception>(comm, error, operation);
        if (global_flags[1] == 0)
          break;
        // ExchangeRecords throws std::logic_error on every rank together, so raising it as
        // `Exception` keeps the ranks in step.
        ReceivedRecords received;
        try
        {
          received = ExchangeRecords(comm, std::move(batch), operation);
        }
        catch (const std::logic_error& exception)
        {
          throw Exception(exception.what());
        }
        try
        {
          source.Receive(stream, pass, received);
        }
        catch (const std::exception& exception)
        {
          error = DescribeError(exception);
        }
      }
    }
  }
}

} // namespace opensn
