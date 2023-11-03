#pragma once

#include <mpi.h>
#include <set>

class Chi;

namespace chi
{

//################################################################### Class def
/**An object for storing various MPI states.*/
class MPI_Info
{
  friend class ::Chi;

private:
  MPI_Comm communicator_ = MPI_COMM_WORLD;
  int location_id_ = 0;
  int process_count_ = 1;

  bool location_id_set_ = false;
  bool process_count_set_ = false;

public:
  const int& location_id = location_id_;     ///< Current process rank.
  const int& process_count = process_count_; ///< Total number of processes.
  const MPI_Comm& comm = communicator_;      ///< MPI communicator

private:
  MPI_Info() = default;

public:
  static MPI_Info& GetInstance() noexcept;

public:
  MPI_Info(const MPI_Info&) = delete;           // Deleted copy constructor
  MPI_Info operator=(const MPI_Info&) = delete; // Deleted assigment operator

protected:
  /**Sets the active communicator*/
  void SetCommunicator(MPI_Comm new_communicator);
  /**Sets the rank.*/
  void SetLocationID(int in_location_id);
  /**Sets the number of processes in the communicator.*/
  void SetProcessCount(int in_process_count);

public:
  /**Calls the generic `MPI_Barrier` with the current communicator.*/
  void Barrier() const;
  static void Call(int mpi_error_code);
};

} // namespace chi
