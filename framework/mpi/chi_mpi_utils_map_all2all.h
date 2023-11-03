#pragma once

#include <map>
#include <vector>

#include <type_traits>

#include "chi_runtime.h"
#include "chi_mpi.h"

namespace chi_mpi_utils
{

/**Given a map with keys indicating the destination process-ids and the
 * values for each key a list of values of type T (T must have an MPI_Datatype).
 * Returns a map with the keys indicating the source process-ids and the
 * values for each key a list of values of type T (sent by the respective
 * process).
 *
 * The keys must be "castable" to `int`.
 *
 * Also expects the MPI_Datatype of T.*/
template <typename K, class T>
std::map<K, std::vector<T>>
MapAllToAll(const std::map<K, std::vector<T>>& pid_data_pairs,
            const MPI_Datatype data_mpi_type,
            const MPI_Comm communicator = Chi::mpi.comm)
{
  static_assert(std::is_integral<K>::value, "Integral datatype required.");

  //============================================= Make sendcounts and
  //                                              senddispls
  std::vector<int> sendcounts(Chi::mpi.process_count, 0);
  std::vector<int> senddispls(Chi::mpi.process_count, 0);
  {
    size_t accumulated_displ = 0;
    for (const auto& [pid, data] : pid_data_pairs)
    {
      sendcounts[pid] = static_cast<int>(data.size());
      senddispls[pid] = static_cast<int>(accumulated_displ);
      accumulated_displ += data.size();
    }
  }

  //============================================= Communicate sendcounts to
  //                                              get recvcounts
  std::vector<int> recvcounts(Chi::mpi.process_count, 0);

  MPI_Alltoall(sendcounts.data(), // sendbuf
               1,
               MPI_INT,           // sendcount, sendtype
               recvcounts.data(), // recvbuf
               1,
               MPI_INT,       // recvcount, recvtype
               communicator); // communicator

  //============================================= Populate recvdispls,
  //                                              sender_pids_set, and
  //                                              total_recv_count
  // All three these quantities are constructed
  // from recvcounts.
  std::vector<int> recvdispls(Chi::mpi.process_count, 0);
  std::set<K> sender_pids_set; // set of neighbor-partitions sending data
  size_t total_recv_count;
  {
    int displacement = 0;
    for (int pid = 0; pid < Chi::mpi.process_count; ++pid)
    {
      recvdispls[pid] = displacement;
      displacement += recvcounts[pid];

      if (recvcounts[pid] > 0) sender_pids_set.insert(static_cast<K>(pid));
    } // for pid
    total_recv_count = displacement;
  }

  //============================================= Make sendbuf
  // The data for each partition is now loaded
  // into a single buffer
  std::vector<T> sendbuf;
  for (const auto& pid_data_pair : pid_data_pairs)
    sendbuf.insert(sendbuf.end(), pid_data_pair.second.begin(), pid_data_pair.second.end());

  //============================================= Make recvbuf
  std::vector<T> recvbuf(total_recv_count);

  //============================================= Communicate serial data
  MPI_Alltoallv(sendbuf.data(),    // sendbuf
                sendcounts.data(), // sendcounts
                senddispls.data(), // senddispls
                data_mpi_type,     // sendtype
                recvbuf.data(),    // recvbuf
                recvcounts.data(), // recvcounts
                recvdispls.data(), // recvdispls
                data_mpi_type,     // recvtype
                communicator);     // comm

  std::map<K, std::vector<T>> output_data;
  {
    for (K pid : sender_pids_set)
    {
      const int data_count = recvcounts.at(pid);
      const int data_displ = recvdispls.at(pid);

      auto& data = output_data[pid];
      data.resize(data_count);

      for (int i = 0; i < data_count; ++i)
        data.at(i) = recvbuf.at(data_displ + i);
    }
  }

  return output_data;
}

} // namespace chi_mpi_utils
