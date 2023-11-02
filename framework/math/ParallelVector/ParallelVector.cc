#include "ParallelVector.h"

#include "mpi/chi_mpi_utils.h"

namespace chi_math
{

ParallelVector::ParallelVector(uint64_t local_size,
                               uint64_t global_size,
                               const MPI_Comm communicator)
  : local_size_(local_size),
    global_size_(global_size),
    location_id_(chi_mpi_utils::GetLocationID(communicator)),
    process_count_(chi_mpi_utils::GetProcessCount(communicator)),
    comm_(communicator)
{
}

ParallelVector::ParallelVector(const ParallelVector& other)
  : local_size_(other.local_size_),
    global_size_(other.global_size_),
    location_id_(other.location_id_),
    process_count_(other.process_count_),
    comm_(other.comm_)
{
}

ParallelVector::ParallelVector(ParallelVector&& other) noexcept
  : local_size_(other.local_size_),
    global_size_(other.global_size_),
    location_id_(other.location_id_),
    process_count_(other.process_count_),
    comm_(other.comm_)
{
}

} // namespace chi_math