// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mpi/record_exchange.h"
#include <algorithm>
#include <limits>
#include <utility>

namespace opensn
{

namespace
{

std::vector<int>
MakeOffsets(const std::vector<int>& counts)
{
  std::vector<int> offsets(counts.size(), 0);
  for (size_t p = 1; p < counts.size(); ++p)
    offsets[p] = offsets[p - 1] + counts[p - 1];
  return offsets;
}

} // namespace

std::int64_t*
RecordBatch::AddRecord(const int owner)
{
  owners.push_back(owner);
  headers.resize(headers.size() + header_width, 0);
  payload_offsets.push_back(payload.size());
  return headers.data() + headers.size() - header_width;
}

void
RecordBatch::AppendPayload(const double* values, const std::size_t count)
{
  payload.insert(payload.end(), values, values + count);
  payload_offsets.back() = payload.size();
}

std::uint64_t
RecordBatch::GetPackedBytes() const
{
  return sizeof(std::int64_t) * (header_width + 1) * GetNumRecords() +
         sizeof(double) * payload.size();
}

std::uint64_t
RecordBatch::GetRecordBytes(const std::uint64_t payload_size) const
{
  return sizeof(std::int64_t) * (header_width + 1) + sizeof(double) * payload_size;
}

ReceivedRecords
ExchangeRecords(const mpi::Communicator& comm, RecordBatch&& records, const std::string& operation)
{
  // Take the records, so that the caller's batch is released as soon as they are packed.
  RecordBatch batch = std::move(records);
  const auto num_ranks = static_cast<size_t>(comm.size());
  const auto num_records = batch.GetNumRecords();
  // Each record carries its payload length after its header values.
  const auto packed_width = batch.header_width + 1;

  // Count header and payload values per destination rank.
  std::string error;
  std::vector<std::uint64_t> send_counts(2 * num_ranks, 0);
  for (size_t r = 0; r < num_records; ++r)
  {
    if (batch.owners[r] < 0 or std::cmp_greater_equal(batch.owners[r], num_ranks))
    {
      error = "a record is addressed to rank " + std::to_string(batch.owners[r]) +
              ", which is not a rank";
      break;
    }
    const auto owner = static_cast<size_t>(batch.owners[r]);
    send_counts[2 * owner] += packed_width;
    send_counts[2 * owner + 1] += batch.payload_offsets[r + 1] - batch.payload_offsets[r];
  }
  ThrowIfAnyRankFailed<std::logic_error>(comm, error, operation);
  std::vector<std::uint64_t> recv_counts(2 * num_ranks, 0);
  comm.all_to_all(send_counts.data(), 2, recv_counts.data(), 2);

  // A rank receives from every rank, so check the totals collectively before using int counts.
  std::array<std::uint64_t, 4> totals = {0, 0, 0, 0};
  for (size_t p = 0; p < num_ranks; ++p)
  {
    totals[0] += send_counts[2 * p];
    totals[1] += send_counts[2 * p + 1];
    totals[2] += recv_counts[2 * p];
    totals[3] += recv_counts[2 * p + 1];
  }
  const auto int_max = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  ThrowIfAnyRankFailed<std::logic_error>(
    comm,
    *std::max_element(totals.begin(), totals.end()) > int_max
      ? "one exchange exceeds the MPI int count range of a rank"
      : "",
    operation);

  std::vector<int> header_counts(num_ranks), payload_counts(num_ranks);
  std::vector<int> recv_header_counts(num_ranks), recv_payload_counts(num_ranks);
  for (size_t p = 0; p < num_ranks; ++p)
  {
    header_counts[p] = static_cast<int>(send_counts[2 * p]);
    payload_counts[p] = static_cast<int>(send_counts[2 * p + 1]);
    recv_header_counts[p] = static_cast<int>(recv_counts[2 * p]);
    recv_payload_counts[p] = static_cast<int>(recv_counts[2 * p + 1]);
  }
  const auto header_offsets = MakeOffsets(header_counts);
  const auto payload_offsets = MakeOffsets(payload_counts);

  // Pack in destination order, keeping the send order of each destination.
  std::vector<std::int64_t> headers(totals[0]);
  std::vector<double> payload(totals[1]);
  {
    auto header_cursor = header_offsets;
    auto payload_cursor = payload_offsets;
    for (size_t r = 0; r < num_records; ++r)
    {
      const auto owner = static_cast<size_t>(batch.owners[r]);
      const auto begin = batch.payload_offsets[r];
      const auto size = batch.payload_offsets[r + 1] - begin;
      // Pointer arithmetic rather than indexing, since a record may have no header values or no
      // payload and the vectors may then be empty.
      auto* h = headers.data() + header_cursor[owner];
      std::copy_n(batch.headers.data() + r * batch.header_width, batch.header_width, h);
      h[batch.header_width] = static_cast<std::int64_t>(size);
      header_cursor[owner] += static_cast<int>(packed_width);
      std::copy_n(batch.payload.data() + begin, size, payload.data() + payload_cursor[owner]);
      payload_cursor[owner] += static_cast<int>(size);
    }
  }
  const auto header_width = batch.header_width;
  batch = RecordBatch(header_width);

  std::vector<std::int64_t> recv_headers;
  comm.all_to_all(headers,
                  header_counts,
                  header_offsets,
                  recv_headers,
                  recv_header_counts,
                  MakeOffsets(recv_header_counts));
  headers = {};
  ReceivedRecords received;
  received.header_width = header_width;
  comm.all_to_all(payload,
                  payload_counts,
                  payload_offsets,
                  received.payload,
                  recv_payload_counts,
                  MakeOffsets(recv_payload_counts));
  payload = {};

  // Unpack in arrival order: by source rank, then in send order.
  const auto num_received = recv_headers.size() / packed_width;
  received.sources.reserve(num_received);
  received.headers.reserve(num_received * header_width);
  received.payload_offsets.reserve(num_received + 1);
  size_t h = 0;
  for (size_t p = 0; p < num_ranks and error.empty(); ++p)
  {
    const auto end = h + static_cast<size_t>(recv_header_counts[p]);
    for (; h < end; h += packed_width)
    {
      const auto size = recv_headers[h + header_width];
      if (size < 0 or received.payload_offsets.back() + static_cast<std::uint64_t>(size) >
                        received.payload.size())
      {
        error = "received record payloads that do not match their headers";
        break;
      }
      received.sources.push_back(static_cast<int>(p));
      received.headers.insert(received.headers.end(),
                              recv_headers.begin() + static_cast<std::ptrdiff_t>(h),
                              recv_headers.begin() + static_cast<std::ptrdiff_t>(h + header_width));
      received.payload_offsets.push_back(received.payload_offsets.back() +
                                         static_cast<std::uint64_t>(size));
    }
  }
  if (error.empty() and received.payload_offsets.back() != received.payload.size())
    error = "received record payloads that do not match their headers";
  ThrowIfAnyRankFailed<std::logic_error>(comm, error, operation);
  return received;
}

OwnerDirectory::OwnerDirectory(const mpi::Communicator& comm,
                               const std::vector<std::uint64_t>& local_keys)
  : comm_(comm)
{
  const auto num_ranks = static_cast<std::uint64_t>(comm.size());
  RecordBatch registrations(1);
  for (const auto key : local_keys)
    registrations.AddRecord(static_cast<int>(key % num_ranks))[0] = static_cast<std::int64_t>(key);
  const auto received =
    ExchangeRecords(comm, std::move(registrations), "Building an owner directory");

  entries_.reserve(received.GetNumRecords());
  for (size_t r = 0; r < received.GetNumRecords(); ++r)
    entries_.emplace_back(static_cast<std::uint64_t>(received.GetHeader(r)[0]),
                          received.sources[r]);
  std::sort(entries_.begin(), entries_.end());
  std::string error;
  for (size_t i = 1; i < entries_.size() and error.empty(); ++i)
    if (entries_[i].first == entries_[i - 1].first)
      error = "key " + std::to_string(entries_[i].first) + " is registered by ranks " +
              std::to_string(entries_[i - 1].second) + " and " + std::to_string(entries_[i].second);
  ThrowIfAnyRankFailed<std::logic_error>(comm, error, "Building an owner directory");
}

std::vector<int>
OwnerDirectory::FindOwners(const std::vector<std::uint64_t>& keys) const
{
  const std::string operation = "Looking up key owners";
  const auto num_ranks = static_cast<std::uint64_t>(comm_.size());

  // Ask the rank that keeps each key, with the key's index so the answer can be placed.
  RecordBatch queries(2);
  for (size_t i = 0; i < keys.size(); ++i)
  {
    auto* header = queries.AddRecord(static_cast<int>(keys[i] % num_ranks));
    header[0] = static_cast<std::int64_t>(keys[i]);
    header[1] = static_cast<std::int64_t>(i);
  }
  const auto asked = ExchangeRecords(comm_, std::move(queries), operation);

  RecordBatch answers(2);
  for (size_t r = 0; r < asked.GetNumRecords(); ++r)
  {
    const auto key = static_cast<std::uint64_t>(asked.GetHeader(r)[0]);
    const auto entry = std::lower_bound(entries_.begin(),
                                        entries_.end(),
                                        key,
                                        [](const auto& e, const auto k) { return e.first < k; });
    auto* header = answers.AddRecord(asked.sources[r]);
    header[0] = asked.GetHeader(r)[1];
    header[1] = entry != entries_.end() and entry->first == key ? entry->second : -1;
  }
  const auto answered = ExchangeRecords(comm_, std::move(answers), operation);

  std::vector<int> owners(keys.size(), -1);
  for (size_t r = 0; r < answered.GetNumRecords(); ++r)
    owners[static_cast<size_t>(answered.GetHeader(r)[0])] =
      static_cast<int>(answered.GetHeader(r)[1]);
  return owners;
}

} // namespace opensn
