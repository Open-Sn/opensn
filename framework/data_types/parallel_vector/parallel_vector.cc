// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/data_types/parallel_vector/parallel_vector.h"

namespace opensn
{

ParallelVector::ParallelVector(uint64_t local_size,
                               uint64_t global_size,
                               const mpi::Communicator& communicator)
  : local_size_(local_size),
    global_size_(global_size),
    location_id_(communicator.rank()),
    process_count_(communicator.size()),
    comm_(communicator)
{
}

ParallelVector::ParallelVector(const ParallelVector& other) = default;

ParallelVector::ParallelVector(ParallelVector&& other) noexcept
  : local_size_(other.local_size_),
    global_size_(other.global_size_),
    location_id_(other.location_id_),
    process_count_(other.process_count_),
    comm_(other.comm_)
{
}

} // namespace opensn
