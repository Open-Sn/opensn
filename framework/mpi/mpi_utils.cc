// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mpi/mpi_utils.h"

namespace opensn
{

std::vector<uint64_t>
BuildLocationExtents(uint64_t local_size, const mpi::Communicator& comm)
{
  const int process_count = comm.size();
  // Get the local vector sizes per process
  std::vector<uint64_t> local_sizes;
  comm.all_gather(local_size, local_sizes);

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  std::vector<uint64_t> extents(process_count + 1, 0);
  for (size_t locJ = 1; locJ < process_count; ++locJ)
    extents[locJ] = extents[locJ - 1] + local_sizes[locJ - 1];
  extents[process_count] = extents[process_count - 1] + local_sizes.back();

  return extents;
}

} // namespace opensn
