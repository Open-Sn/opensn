#include "framework/mpi/mpi_info.h"

#include "framework/logging/chi_log_exceptions.h"

namespace chi
{

MPI_Info&
MPI_Info::GetInstance() noexcept
{
  static MPI_Info singleton;
  return singleton;
}

void
MPI_Info::SetCommunicator(MPI_Comm new_communicator)
{
  communicator_ = new_communicator;
}

void
MPI_Info::SetLocationID(int in_location_id)
{
  if (not location_id_set_) location_id_ = in_location_id;
  location_id_set_ = true;
}

void
MPI_Info::SetProcessCount(int in_process_count)
{
  if (not process_count_set_) process_count_ = in_process_count;
  process_count_set_ = true;
}

void
MPI_Info::Barrier() const
{
  MPI_Barrier(this->communicator_);
}

void
MPI_Info::Call(int mpi_error_code)
{
  if (mpi_error_code == MPI_SUCCESS) return;
  else
  {
    char estring[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string(mpi_error_code, estring, &resultlen);

    ChiLogicalError(estring);
  }
}

} // namespace chi
