#include "framework/math/ParallelVector/ParallelSTLVector.h"

#include "framework/mpi/chi_mpi_utils_map_all2all.h"
#include "framework/data_types/byte_array.h"

#include <petsc.h>

#include "framework/logging/chi_log.h"
#include "framework/logging/chi_log_exceptions.h"

#include <sstream>
#include <stdexcept>
#include <cmath>
#include <numeric>

#define scint64_t static_cast<int64_t>
#define scint static_cast<int>
#define scshort static_cast<short>

namespace chi_math
{

ParallelSTLVector::ParallelSTLVector(const uint64_t local_size,
                                     const uint64_t global_size,
                                     const MPI_Comm communicator)
  : ParallelVector(local_size, global_size, communicator),
    extents_(DefineExtents(local_size, process_count_, communicator)),
    values_(local_size_, 0.0)
{
}

ParallelSTLVector::ParallelSTLVector(const ParallelSTLVector& other)
  : ParallelVector(other), extents_(other.extents_), values_(other.values_)
{
}

ParallelSTLVector::ParallelSTLVector(ParallelSTLVector&& other) noexcept
  : ParallelVector(std::move(other)), extents_(other.extents_), values_(std::move(other.values_))
{
}

std::unique_ptr<ParallelVector>
ParallelSTLVector::MakeCopy() const
{
  return std::make_unique<ParallelSTLVector>(*this);
}

std::unique_ptr<ParallelVector>
ParallelSTLVector::MakeClone() const
{
  auto new_vec = std::make_unique<ParallelSTLVector>(*this);

  new_vec->Set(0.0);

  return new_vec;
}

double*
ParallelSTLVector::Data()
{
  return values_.data();
}
const double*
ParallelSTLVector::Data() const
{
  return values_.data();
}

const std::vector<double>&
ParallelSTLVector::LocalSTLData() const
{
  return values_;
}

std::vector<double>&
ParallelSTLVector::LocalSTLData()
{
  return values_;
}

std::vector<double>
ParallelSTLVector::MakeLocalVector()
{
  return std::vector<double>(values_.begin(), values_.begin() + scint(local_size_));
}

double
ParallelSTLVector::operator[](const int64_t local_id) const
{
  ChiInvalidArgumentIf(local_id < 0 or local_id >= values_.size(),
                       "Invalid local index provided. " + std::to_string(local_id) + " vs [0," +
                         std::to_string(values_.size()) + ")");

  return values_[local_id];
}

double&
ParallelSTLVector::operator[](const int64_t local_id)
{
  ChiInvalidArgumentIf(local_id < 0 or local_id >= values_.size(),
                       "Invalid local index provided. " + std::to_string(local_id) + " vs [0," +
                         std::to_string(values_.size()) + ")");

  return values_[local_id];
}

void
ParallelSTLVector::Set(const double value)
{
  values_.assign(values_.size(), value);
}

void
ParallelSTLVector::Set(const std::vector<double>& local_vector)
{
  ChiInvalidArgumentIf(local_vector.size() < local_size_, "Incompatible local vector size.");

  // We cannot assign values = local_vector because this might
  // destroy the internals of a ghosted vector
  std::copy(local_vector.begin(), local_vector.begin() + scint64_t(local_size_), values_.begin());
}

void
ParallelSTLVector::BlockSet(const std::vector<double>& y, int64_t local_offset, int64_t num_values)
{
  ChiInvalidArgumentIf(y.size() < num_values,
                       "y.size() < num_values " + std::to_string(y.size()) + " < " +
                         std::to_string(num_values));

  const int64_t local_end = local_offset + num_values;
  ChiInvalidArgumentIf(local_end > local_size_,
                       "local_offset + num_values=" + std::to_string(local_end) +
                         ", is out of range for destination vector with local size " +
                         std::to_string(local_size_));

  std::copy(y.begin(), y.begin() + num_values, values_.begin() + local_offset);
}

void
ParallelSTLVector::CopyLocalValues(const ParallelVector& y)
{
  ChiLogicalErrorIf(y.LocalSize() != local_size_, "y.LocalSize() != local_size_");
  ChiLogicalErrorIf(y.GlobalSize() != global_size_, "y.GlobalSize() != global_size_");

  // This prevents us from unnecessarily going through the operator[]
  // of y, which again applies bounds checking
  const double* y_data = y.Data();

  std::copy(y_data, y_data + local_size_, values_.begin());
}

void
ParallelSTLVector::CopyLocalValues(Vec y)
{
  PetscInt n;
  VecGetLocalSize(y, &n);

  ChiInvalidArgumentIf(n < local_size_, "Attempted update with a vector of insufficient size.");

  const double* x;
  VecGetArrayRead(y, &x);
  std::copy(x, x + n, values_.begin());
  VecRestoreArrayRead(y, &x);
}

void
ParallelSTLVector::BlockCopyLocalValues(const ParallelVector& y,
                                        int64_t y_offset,
                                        int64_t local_offset,
                                        int64_t num_values)
{
  ChiInvalidArgumentIf(y_offset < 0, "y_offset < 0.");
  ChiInvalidArgumentIf(local_offset < 0, "local_offset < 0.");

  const int64_t y_end = y_offset + num_values;
  const int64_t local_end = local_offset + num_values;

  ChiInvalidArgumentIf(y_end > y.LocalSize(),
                       "y_offset + num_values=" + std::to_string(y_end) +
                         ", is out of range for vector y with local size " +
                         std::to_string(y.LocalSize()));

  ChiInvalidArgumentIf(local_end > local_size_,
                       "local_offset + num_values=" + std::to_string(local_end) +
                         ", is out of range for destination vector with local size " +
                         std::to_string(local_size_));

  // This prevents us from unnecessarily going through the operator[]
  // of y, which again applies bounds checking
  const double* y_data = y.Data();

  std::copy(y_data + y_offset, y_data + y_offset + num_values, values_.begin() + local_offset);
}

void
ParallelSTLVector::BlockCopyLocalValues(Vec y,
                                        int64_t y_offset,
                                        int64_t local_offset,
                                        int64_t num_values)
{
  ChiInvalidArgumentIf(y_offset < 0, "y_offset < 0.");
  ChiInvalidArgumentIf(local_offset < 0, "local_offset < 0.");

  const int64_t y_end = y_offset + num_values;
  const int64_t local_end = local_offset + num_values;

  int64_t y_local_size;
  VecGetLocalSize(y, &y_local_size);

  ChiInvalidArgumentIf(y_end > y_local_size,
                       "y_offset + num_values=" + std::to_string(y_end) +
                         ", is out of range for vector y with local size " +
                         std::to_string(y_local_size));

  ChiInvalidArgumentIf(local_end > local_size_,
                       "local_offset + num_values=" + std::to_string(local_end) +
                         ", is out of range for destination vector with local size " +
                         std::to_string(local_size_));

  const double* y_data;
  VecGetArrayRead(y, &y_data);

  std::copy(y_data + y_offset, y_data + y_offset + num_values, values_.begin() + local_offset);

  VecRestoreArrayRead(y, &y_data);
}

void
ParallelSTLVector::SetValue(const int64_t global_id, const double value, const VecOpType op_type)
{
  ChiInvalidArgumentIf(global_id < 0 or global_id >= global_size_,
                       "Invalid global index encountered. Global indices "
                       "must be in the range [0, this->GlobalSize()].");

  auto& op_cache = op_type == VecOpType::SET_VALUE ? set_cache_ : add_cache_;
  op_cache.emplace_back(global_id, value);
}

void
ParallelSTLVector::SetValues(const std::vector<int64_t>& global_ids,
                             const std::vector<double>& values,
                             const VecOpType op_type)
{
  ChiInvalidArgumentIf(global_ids.size() != values.size(),
                       "Size mismatch between indices and values.");

  auto& op_cache = op_type == VecOpType::SET_VALUE ? set_cache_ : add_cache_;
  for (size_t i = 0; i < global_ids.size(); ++i)
  {
    const auto& global_id = global_ids[i];
    ChiInvalidArgumentIf(global_id < 0 or global_id >= global_size_,
                         "Invalid global index encountered. Global indices "
                         "must be in the range [0, this->GlobalSize()].");
    op_cache.emplace_back(global_id, values[i]);
  }
}

void
ParallelSTLVector::operator+=(const ParallelVector& y)
{
  ChiLogicalErrorIf(y.LocalSize() != local_size_, "y.LocalSize() != local_size_");
  ChiLogicalErrorIf(y.GlobalSize() != global_size_, "y.GlobalSize() != global_size_");

  // This prevents us from unnecessarily going through the operator[]
  // of y, which again applies bounds checking
  const double* y_data = y.Data();

  for (int64_t i = 0; i < local_size_; ++i)
    values_[i] += y_data[i];
}

void
ParallelSTLVector::PlusAY(const ParallelVector& y, double a)
{
  ChiLogicalErrorIf(y.LocalSize() != local_size_, "y.LocalSize() != local_size_");
  ChiLogicalErrorIf(y.GlobalSize() != global_size_, "y.GlobalSize() != global_size_");

  // This prevents us from unnecessarily going through the operator[]
  // of y, which again applies bounds checking
  const double* y_data = y.Data();

  if (a == 1.0)
    for (int64_t i = 0; i < local_size_; ++i)
      values_[i] += y_data[i];
  else if (a == -1.0)
    for (int64_t i = 0; i < local_size_; ++i)
      values_[i] -= y_data[i];
  else
    for (int64_t i = 0; i < local_size_; ++i)
      values_[i] += a * y_data[i];
}

void
ParallelSTLVector::AXPlusY(double a, const ParallelVector& y)
{
  ChiLogicalErrorIf(y.LocalSize() != local_size_, "y.LocalSize() != local_size_");
  ChiLogicalErrorIf(y.GlobalSize() != global_size_, "y.GlobalSize() != global_size_");

  // This prevents us from unnecessarily going through the operator[]
  // of y, which again applies bounds checking
  const double* y_data = y.Data();

  for (int64_t i = 0; i < local_size_; ++i)
    values_[i] = a * values_[i] + y_data[i];
}

void
ParallelSTLVector::Scale(double a)
{
  for (size_t i = 0; i < local_size_; ++i)
    values_[i] *= a;
}

void
ParallelSTLVector::Shift(double a)
{
  for (size_t i = 0; i < local_size_; ++i)
    values_[i] += a;
}

/**Returns the specified norm of the vector.*/
double
ParallelSTLVector::ComputeNorm(chi_math::NormType norm_type) const
{
  switch (norm_type)
  {
    case NormType::L1_NORM:
    {
      const double local_norm_val =
        std::accumulate(values_.begin(), values_.begin() + scint(local_size_), 0.0);

      double global_norm_val;
      MPI_Allreduce(&local_norm_val,  // sendbuf
                    &global_norm_val, // recvbuf
                    1,                // sendcount
                    MPI_DOUBLE,       // sendtype
                    MPI_SUM,          // operation
                    comm_);           // communicator

      return global_norm_val;
    }
    case NormType::L2_NORM:
    {
      double local_norm_val = 0.0;
      for (size_t i = 0; i < local_size_; ++i)
      {
        const double value = values_[i];
        local_norm_val += value * value;
      }
      double global_norm_val;
      MPI_Allreduce(&local_norm_val,  // sendbuf
                    &global_norm_val, // recvbuf
                    1,                // sendcount
                    MPI_DOUBLE,       // sendtype
                    MPI_SUM,          // operation
                    comm_);           // communicator

      return std::sqrt(global_norm_val);
    }
    case NormType::LINF_NORM:
    {
      const double local_norm_val =
        *std::max_element(values_.begin(), values_.begin() + scint(local_size_));

      double global_norm_val;
      MPI_Allreduce(&local_norm_val,  // sendbuf
                    &global_norm_val, // recvbuf
                    1,                // sendcount
                    MPI_DOUBLE,       // sendtype
                    MPI_MAX,          // operation
                    comm_);           // communicator

      return global_norm_val;
    }
    default:
      return 0.0;
  }
}

void
ParallelSTLVector::Assemble()
{
  // Define the local operation mode.
  // 0=Do Nothing, 1=Set, 2=Add, 3=INVALID (mixed set/add ops)
  const short local_mode =
    scshort(not set_cache_.empty()) + scshort(not add_cache_.empty()) * short(2);
  ChiLogicalErrorIf(local_mode == 3, "Invalid operation mode.");

  // Now, determine the global operation mode
  short global_mode;
  MPI_Allreduce(&local_mode, &global_mode, 1, MPI_SHORT, MPI_MAX, comm_);

  // If the mode is to do nothing, exit
  if (global_mode == 0) return;

  // Next, ensure that all operation types are compatible
  ChiLogicalErrorIf(local_mode != 0 and local_mode != global_mode,
                    "The operation on each process must be either 0 (do nothing),"
                    "or the same across all processes.");

  // Now, store the global operation type and get the appropriate cache
  using OpType = VecOpType;
  const auto global_op_type = static_cast<OpType>(global_mode);
  auto& op_cache = global_op_type == OpType ::SET_VALUE ? set_cache_ : add_cache_;

  // First, segregate the local and non-local operations
  std::vector<std::pair<int64_t, double>> local_cache;
  std::vector<std::pair<int64_t, double>> nonlocal_cache;
  for (const auto& op : op_cache)
  {
    const int op_pid = FindOwnerPID(op.first);
    if (op_pid == location_id_) local_cache.emplace_back(op);
    else
      nonlocal_cache.emplace_back(op);
  }

  // The local operations can be handled immediately
  for (const auto& [global_id, value] : local_cache)
  {
    const int64_t local_id = global_id - scint64_t(extents_[location_id_]);
    ChiLogicalErrorIf(local_id < 0 or local_id >= local_size_,
                      "Invalid mapping from global to local.");

    if (global_op_type == OpType::SET_VALUE) values_[local_id] = value;
    else
      values_[local_id] += value;
  }

  // With this, the data that needs to be sent to other processes must be
  // determined. Here, a mapping is developed between the processes that
  // need to be sent information, and the serialized operations that need
  // to be sent. The operations are serialized by converting the
  // int64_t-double pair to bytes.
  std::map<int, chi_data_types::ByteArray> pid_send_map;
  for (const auto& [global_id, value] : nonlocal_cache)
  {
    const int pid = FindOwnerPID(global_id);
    auto& byte_array = pid_send_map[pid];
    byte_array.Write(global_id);
    byte_array.Write(value);
  }

  // For use with MPI, the byte arrays from above is converted to an
  // STL vector of bytes.
  std::map<int, std::vector<std::byte>> pid_send_map_bytes;
  for (const auto& [pid, byte_array] : pid_send_map)
    pid_send_map_bytes[pid] = byte_array.Data();

  // With the information that needs to be sent to other processes obtained,
  // now, the information to be received from other processes is needed.
  // To do this, each process must send to each other process the information
  // that it needs. With each process knowing what each other process needs
  // from it, a map of information to be sent is obtained.
  std::map<int, std::vector<std::byte>> pid_recv_map_bytes =
    chi_mpi_utils::MapAllToAll(pid_send_map_bytes, MPI_BYTE);

  // The received information is now processed, unpacked, and the
  // necessary operations performed
  for (const auto& [pid, byte_vector] : pid_recv_map_bytes)
  {
    const auto packet_size = sizeof(std::pair<int64_t, double>);

    ChiLogicalErrorIf(byte_vector.size() % packet_size != 0,
                      "Unrecognized received operations. Operations are serialized with "
                      "an int64_t and double, but the received packet from process " +
                        std::to_string(pid) + " on process " + std::to_string(location_id_) +
                        " is not an integer multiple of the size of an int64_t and double.");

    const size_t num_ops = byte_vector.size() / packet_size;
    chi_data_types::ByteArray byte_array(byte_vector);
    for (size_t k = 0; k < num_ops; ++k)
    {
      const int64_t global_id = byte_array.Read<int64_t>();
      const double value = byte_array.Read<double>();

      // Check that the global ID is in fact valid for this process
      const int64_t local_id = global_id - scint64_t(extents_[location_id_]);

      ChiLogicalErrorIf(local_id < 0 or local_id >= local_size_,
                        "A non-local global ID was received by process " +
                          std::to_string(location_id_) + " by process " + std::to_string(pid) +
                          " during vector assembly.");

      // Contribute to the local vector
      if (global_op_type == OpType ::SET_VALUE) values_[local_id] = value;
      else
        values_[local_id] += value;
    }
  }

  // Finally, clear the operation cache
  op_cache.clear();
}

std::vector<uint64_t>
ParallelSTLVector::DefineExtents(uint64_t local_size, int comm_size, MPI_Comm communicator)
{
  // Get the local vector sizes per processor
  std::vector<uint64_t> local_sizes(comm_size, 0);
  MPI_Allgather(&local_size,        // sendbuf
                1,                  // sendcount
                MPI_UINT64_T,       // sendcount
                local_sizes.data(), // recvbuf
                1,                  // recvcount
                MPI_UINT64_T,       // recvtype
                communicator);      // communicator

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  std::vector<uint64_t> extents(comm_size + 1, 0);
  for (size_t p = 1; p < comm_size; ++p)
    extents[p] = extents[p - 1] + local_sizes[p - 1];
  extents[comm_size] = extents[comm_size - 1] + local_sizes.back();

  return extents;
}

int
ParallelSTLVector::FindOwnerPID(const uint64_t global_id) const
{
  ChiInvalidArgumentIf(global_id >= global_size_, "Invalid global id specified.");

  for (int p = 0; p < process_count_; ++p)
    if (global_id >= extents_[p] and global_id < extents_[p + 1]) return p;
  return -1;
}

std::string
ParallelSTLVector::PrintStr() const
{
  std::stringstream ss;
  for (size_t i = 0; i < values_.size(); ++i)
    ss << (i == 0 ? "[" : "") << values_[i] << (i < values_.size() - 1 ? " " : "]");
  return ss.str();
}

} // namespace chi_math
