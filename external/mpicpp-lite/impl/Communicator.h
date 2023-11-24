#pragma once

#include "mpi.h"
#include <vector>
#include <cassert>
#include "Datatype.h"
#include "Status.h"
#include "Request.h"
#include "Operation.h"
#include "Error.h"
#include "Group.h"

namespace mpicpp_lite {

/// Wrapper around `MPI_Comm`
class Communicator {
public:
    /// Create `MPI_COMM_WORLD` communicator
    Communicator();

    /// Create communicator from an `MPI_Comm` one
    Communicator(const MPI_Comm & comm);

    /// Copy constructor
    Communicator(const Communicator & comm);

    /// Determine the rank of the executing process in a communicator
    ///
    /// @return Rank of the executing process
    int rank() const;

    /// Determine the number of processes in a communicator
    ///
    /// @return Number of processes
    int size() const;

    ///
    Communicator create(const Group & group, int tag = 0) const;

    /// Accesses the group associated with given communicator
    ///
    /// @return Group corresponding to communicator
    Group group() const;

    ///
    bool is_valid() const;

    /// Send data to another process
    ///
    /// @tparam T C++ type of the data
    /// @param dest Destination rank
    /// @param tag Message tag
    /// @param value Value to send
    template <typename T>
    void send(int dest, int tag, const T & value) const;

    /// Send data to another process
    ///
    /// @tparam T C++ type of the data
    /// @param dest Destination rank
    /// @param tag Message tag
    /// @param values Values to send
    /// @param n Number of values to send
    template <typename T>
    void send(int dest, int tag, const T * values, int n) const;

    /// Send `std::vector` of data to another process
    ///
    /// @tparam T C++ type of the data
    /// @param dest Destination rank
    /// @param tag Message tag
    /// @param value Vector of `T` to send
    template <typename T, typename A>
    void send(int dest, int tag, const std::vector<T, A> & value) const;

    /// Send a message to another process without any data
    ///
    /// @param dest Destination rank
    /// @param tag Message tag
    void send(int dest, int tag) const;

    /// Receive data from a remote process
    ///
    /// @tparam T C++ type of the data
    /// @param source Source rank
    /// @param tag Message tag
    /// @param value Variable to recieve the data
    /// @return `Status` of the operation
    template <typename T>
    Status recv(int source, int tag, T & value) const;

    /// Receive data from a remote process
    ///
    /// @tparam T C++ type of the data
    /// @param source Source rank
    /// @param tag Message tag
    /// @param values Variable to recieve the data
    /// @param n Number of values to receive
    /// @return `Status` of the operation
    template <typename T>
    Status recv(int source, int tag, T * values, int n) const;

    /// Receive std::vector of data from a remote process
    ///
    /// @tparam T C++ type of the data
    /// @param source Source rank
    /// @param tag Message tag
    /// @param value Variable to recieve the data
    /// @return `Status` of the operation
    template <typename T, typename A>
    Status recv(int source, int tag, std::vector<T, A> & value) const;

    /// Receive a message from a remote process without any data
    ///
    /// @param source Source rank
    /// @param tag Message tag
    /// @return `Status` of the operation
    Status recv(int source, int tag) const;

    /// Send a message to a remote process without blocking
    ///
    /// @tparam T C++ type of the data
    /// @param dest Destination rank
    /// @param tag Message tag
    /// @param value Value to send
    /// @return Communication `Request`
    template <typename T>
    Request isend(int dest, int tag, const T & value) const;

    /// Send a std::vector of values to a remote process without blocking
    ///
    /// @tparam T C++ type of the data
    /// @param dest Destination rank
    /// @param tag Message tag
    /// @param value Values to send
    /// @return Communication `Request`
    template <typename T, typename A>
    Request isend(int dest, int tag, const std::vector<T, A> & values) const;

    /// Send a message to a remote process without blocking
    ///
    /// @tparam T C++ type of the data
    /// @param dest Destination rank
    /// @param tag Message tag
    /// @param values Values to send
    /// @param n Number of values to send
    /// @return Communication `Request`
    template <typename T>
    Request isend(int dest, int tag, const T * values, int n) const;

    /// Receive a message from a remote process without blocking
    ///
    /// @tparam T C++ type of the data
    /// @param source Source rank
    /// @param tag Message tag
    /// @param value Variable to recieve the data
    /// @return Communication `Request`
    template <typename T>
    Request irecv(int source, int tag, T & value) const;

    /// Receive a message from a remote process without blocking
    ///
    /// @tparam T C++ type of the data
    /// @param source Source rank
    /// @param tag Message tag
    /// @param values Variable to recieve the data
    /// @param n Number of values to receive
    template <typename T>
    Request irecv(int source, int tag, T * values, int n) const;

    /// Nonblocking test for a message
    ///
    /// @param source Rank of source or `ANY_SOURCE`
    /// @param tag Message tag or `ANY_TAG`
    /// @return `true` if a message with the specified source, and tag is available, `false`
    ///          otherwise
    bool iprobe(int source, int tag) const;

    /// Nonblocking test for a message
    ///
    /// @param source Rank of source or `ANY_SOURCE`
    /// @param tag Message tag or `ANY_TAG`
    /// @param status Status object
    /// @return `true` if a message with the specified source, and tag is available, `false`
    ///          otherwise
    bool iprobe(int source, int tag, Status & status) const;

    /// Wait for all processes within a communicator to reach the barrier.
    void barrier() const;

    /// Broadcast a value from a root process to all other processes
    ///
    /// @tparam T C++ type of the data
    /// @param value Value to send
    /// @param root Rank of the sending process
    template <typename T>
    void broadcast(T & value, int root) const;

    /// Broadcast a std::vector of values from a root process to all other processes
    ///
    /// @tparam T C++ type of the data
    /// @param value Value to send
    /// @param root Rank of the sending process
    template <typename T>
    void broadcast(std::vector<T> & value, int root) const;

    /// Broadcast a value from a root process to all other processes
    ///
    /// @tparam T C++ type of the data
    /// @param values Values to send
    /// @param n Number of values to send
    /// @param root Rank of the sending process
    template <typename T>
    void broadcast(T * values, int n, int root) const;

    /// Gather together values from a group of processes
    ///
    /// @tparam T C++ type of the data
    /// @param in_value Value to send
    /// @param out_values Receiving variable
    /// @param root Rank of receiving process
    template <typename T>
    void gather(const T & in_value, T * out_values, int root) const;

    /// Gather together values from a group of processes
    ///
    /// @tparam T C++ type of the data
    /// @param in_value Value to send
    /// @param out_values Receiving variable
    /// @param root Rank of receiving process
    template <typename T>
    void gather(const T & in_value, std::vector<T> & out_values, int root) const;

    /// Gather together values from a group of processes
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param n Number of values to send
    /// @param out_values Receiving variable
    /// @param root Rank of receiving process
    template <typename T>
    void gather(const T * in_values, int n, T * out_values, int root) const;

    /// Gather together values from a group of processes
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param n Number of values to send
    /// @param out_values Receiving variable
    /// @param root Rank of receiving process
    template <typename T>
    void gather(const T * in_values, int n, std::vector<T> & out_values, int root) const;

    /// Gathers data from all tasks and distribute the combined data to all tasks
    ///
    /// @tparam T C++ type of the data
    /// @param in_value Value to send
    /// @param n Number of values to send
    /// @param out_values Receiving variable
    /// @param m Number of values to receive
    template <typename T>
    void all_gather(const T * in_value, int n, T * out_values, int m) const;

    /// Gathers data from all tasks and distribute the combined data to all tasks
    ///
    /// @tparam T C++ type of the data
    /// @param in_value Value to send
    /// @param out_values Receiving variable
    template <typename T>
    void all_gather(const T & in_value, std::vector<T> & out_values) const;

    /// Gathers data from all tasks and distribute the combined data to all tasks
    ///
    /// @tparam T C++ type of the data
    /// @param in_value Vector of values to send (the size of this vector can be different on every
    ///                 process)
    /// @param out_values Receiving variable
    template <typename T>
    void all_gather(const std::vector<T> & in_value, std::vector<T> & out_values) const;

    /// Gathers data from all tasks and deliver the combined data to all tasks
    ///
    /// @tparam T C++ datatype
    /// @param in_values[in] Values to send
    /// @param out_values[out] Buffer to receive the data
    /// @param out_counts[in] Integer array (of length group size) containing the number of elements
    ///        that are to be received from each process
    /// @param out_offsets[in] Integer array (of length group size). Entry `i` specifies the
    ///        displacement (relative to out_values) at which to place the incoming data from
    ///        process `i`
    template <typename T>
    void all_gather(const std::vector<T> & in_values,
                    std::vector<T> & out_values,
                    const std::vector<int> & out_counts,
                    const std::vector<int> & out_offsets) const;

    /// Send data from one process to all other processes in a communicator
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param out_value Receiving variable
    /// @param root Rank of the sending process
    template <typename T>
    void scatter(const T * in_values, T & out_value, int root) const;

    /// Send data from one process to all other processes in a communicator
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param out_value Receiving variable
    /// @param root Rank of the sending process
    template <typename T>
    void scatter(const std::vector<T> & in_values, T & out_value, int root) const;

    /// Send data from one process to all other processes in a communicator
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param out_values Receiving variable
    /// @param n Number of values to send
    /// @param root Rank of the sending process
    template <typename T>
    void scatter(const T * in_values, T * out_values, int n, int root) const;

    /// Send data from one process to all other processes in a communicator
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param out_values Receiving variable
    /// @param n Number of values to send
    /// @param root Rank of the sending process
    template <typename T>
    void scatter(const std::vector<T> & in_values, T * out_values, int n, int root) const;

    /// Reduce values on all processes to a single value
    ///
    /// @tparam T C++ type of the data
    /// @tparam Op Type of the reduce operation
    /// @param in_values Values to send
    /// @param n Number of values to send
    /// @param out_values Receiving variable
    /// @param op Reduce operation
    /// @param root Rank of root process
    template <typename T, typename Op>
    void reduce(const T * in_values, int n, T * out_values, Op op, int root) const;

    /// Reduce values on all processes to a single value
    ///
    /// @tparam T C++ type of the data
    /// @tparam Op Type of the reduce operation
    /// @param in_values Values to send
    /// @param out_values Receiving variable
    /// @param op Reduce operation
    /// @param root Rank of root process
    template <typename T, typename Op>
    void
    reduce(std::vector<T> const & in_values, std::vector<T> & out_values, Op op, int root) const;

    /// Reduce values on all processes to a single value
    ///
    /// @tparam T C++ type of the data
    /// @tparam Op Type of the reduce operation
    /// @param in_value Values to send
    /// @param out_value Receiving variable
    /// @param op Reduce operation
    /// @param root Rank of root process
    template <typename T, typename Op>
    void reduce(const T & in_value, T & out_value, Op op, int root) const;

    /// Combine values from all processes and distributes the result back to all processes
    ///
    /// @tparam T C++ type of the data
    /// @tparam Op Type of the reduce operation
    /// @param in_values Values to send
    /// @param n Number of values to send
    /// @param out_values Receiving variable
    /// @param op Reduce operation
    template <typename T, typename Op>
    void all_reduce(const T * in_values, int n, T * out_values, Op op) const;

    /// Combine values from all processes and distributes the result back to all processes
    ///
    /// @tparam T C++ type of the data
    /// @tparam Op Type of the reduce operation
    /// @param in_value Values to send
    /// @param out_value Receiving variable
    /// @param op Reduce operation
    template <typename T, typename Op>
    void all_reduce(const T & in_value, T & out_value, Op op) const;

    /// Combine values from all processes and distributes the result back to all processes
    ///
    /// @tparam T C++ type of the data
    /// @tparam Op Type of the reduce operation
    /// @param value Send/receive variable
    /// @param op Reduce operation
    template <typename T, typename Op>
    void all_reduce(T & value, Op op) const;

    /// Sends data from all to all processes
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param n Number valus to send
    /// @param out_values Receiving variable
    /// @param m Values to receive
    template <typename T>
    void all_to_all(const T * in_values, int n, T * out_values, int m) const;

    /// Sends data from all to all processes
    ///
    /// @tparam T C++ type of the data
    /// @param in_values Values to send
    /// @param out_values Receiving variable
    template <typename T>
    void all_to_all(const std::vector<T> & in_values, std::vector<T> & out_values) const;

    /// Sends data from all to all processes
    template <typename T>
    void all_to_all(const std::vector<std::vector<T>> & in_values,
                    std::vector<T> & out_values) const;

    /// Sends data from all to all processes; each process may send a different amount of data and
    /// provide displacements for the input and output data.
    ///
    /// @tparam T C++ datatype
    /// @param in_values[in] Values to send
    /// @param in_counts[in] Integer array equal to the group size specifying the number of elements
    ///        to send to each processor
    /// @param in_offsets[in] integer array (of length group size). Entry `j` specifies the
    ///        displacement relative to sendbuf from which to take the outgoing data destined for
    ///        process `j`
    /// @param out_values[out] Buffer that will receive the data
    /// @param out_counts[in] Integer array equal to the group size specifying the maximum number of
    ///        elements that can be received from each processor
    /// @param out_offsets[in] Integer array (of length group size). Entry `i` specifies the
    ///        displacement relative to recvbuf at which to place the incoming data from process `i`
    template <typename T>
    void all_to_all(const std::vector<T> & in_values,
                    const std::vector<int> & in_counts,
                    const std::vector<int> & in_offsets,
                    std::vector<T> & out_values,
                    const std::vector<int> & out_counts,
                    const std::vector<int> & out_offsets) const;

    /// Abort all tasks in the group of this communicator
    ///
    /// @param errcode Error code to return to invoking environment
    void abort(int errcode) const;

    /// Cast operator so we can pass this directly in MPI API
    operator MPI_Comm() const;

    void set_error_handler();

private:
    MPI_Comm comm;
};

//

inline Communicator::Communicator() : comm(MPI_COMM_WORLD) {}

inline Communicator::Communicator(const MPI_Comm & comm) : comm(comm) {}

inline Communicator::Communicator(const Communicator & comm) : comm(comm.comm) {}

inline int
Communicator::rank() const
{
    int r;
    MPI_Comm_rank(this->comm, &r);
    return r;
}

inline int
Communicator::size() const
{
    int sz;
    MPI_Comm_size(this->comm, &sz);
    return sz;
}

inline Communicator
Communicator::create(const Group & group, int tag) const
{
    MPI_Comm new_comm;
    MPI_CHECK_SELF(MPI_Comm_create_group(this->comm, group, tag, &new_comm));
    return { new_comm };
}

inline Group
Communicator::group() const
{
    MPI_Group g;
    MPI_CHECK_SELF(MPI_Comm_group(this->comm, &g));
    return { g };
}

inline bool
Communicator::is_valid() const
{
    return this->comm != MPI_COMM_NULL;
}

// Send

template <typename T>
void
Communicator::send(int dest, int tag, const T & value) const
{
    send(dest, tag, &value, 1);
}

template <typename T>
void
Communicator::send(int dest, int tag, const T * values, int n) const
{
    assert(values != nullptr);
    MPI_CHECK_SELF(
        MPI_Send(const_cast<T *>(values), n, get_mpi_datatype<T>(), dest, tag, this->comm));
}

template <typename T, typename A>
void
Communicator::send(int dest, int tag, const std::vector<T, A> & value) const
{
    typename std::vector<T, A>::size_type size = value.size();
    send(dest, tag, value.data(), size);
}

inline void
Communicator::send(int dest, int tag) const
{
    MPI_CHECK_SELF(MPI_Send(MPI_BOTTOM, 0, MPI_PACKED, dest, tag, this->comm));
}

// Recv

template <typename T>
Status
Communicator::recv(int source, int tag, T & value) const
{
    return recv(source, tag, &value, 1);
}

template <typename T>
Status
Communicator::recv(int source, int tag, T * values, int n) const
{
    assert(values != nullptr);
    MPI_Status status = { 0 };
    MPI_CHECK_SELF(MPI_Recv(const_cast<T *>(values),
                            n,
                            get_mpi_datatype<T>(),
                            source,
                            tag,
                            this->comm,
                            &status));
    return { status };
}

template <typename T, typename A>
Status
Communicator::recv(int source, int tag, std::vector<T, A> & values) const
{
    MPI_Status status = { 0 };
    MPI_CHECK_SELF(MPI_Probe(source, tag, this->comm, &status));
    int size = 0;
    MPI_Get_count(&status, get_mpi_datatype<T>(), &size);
    values.resize(size);
    return recv(source, tag, values.data(), size);
}

inline Status
Communicator::recv(int source, int tag) const
{
    MPI_Status status = { 0 };
    MPI_CHECK_SELF(MPI_Recv(MPI_BOTTOM, 0, MPI_PACKED, source, tag, this->comm, &status));
    return { status };
}

// Isend

template <typename T>
Request
Communicator::isend(int dest, int tag, const T & value) const
{
    return isend(dest, tag, &value, 1);
}

template <typename T, typename A>
Request
Communicator::isend(int dest, int tag, const std::vector<T, A> & values) const
{
    return isend(dest, tag, values.data(), values.size());
}

template <typename T>
Request
Communicator::isend(int dest, int tag, const T * values, int n) const
{
    assert(values != nullptr);
    MPI_Request request;
    MPI_CHECK_SELF(MPI_Isend(const_cast<T *>(values),
                             n,
                             get_mpi_datatype<T>(),
                             dest,
                             tag,
                             this->comm,
                             &request));
    return { request };
}

// Irecv

template <typename T>
Request
Communicator::irecv(int source, int tag, T & value) const
{
    return irecv(source, tag, &value, 1);
}

template <typename T>
Request
Communicator::irecv(int source, int tag, T * values, int n) const
{
    assert(values != nullptr);
    MPI_Request request;
    MPI_CHECK_SELF(MPI_Irecv(const_cast<T *>(values),
                             n,
                             get_mpi_datatype<T>(),
                             source,
                             tag,
                             this->comm,
                             &request));
    return { request };
}

inline bool
Communicator::iprobe(int source, int tag) const
{
    int flag;
    MPI_CHECK_SELF(MPI_Iprobe(source, tag, this->comm, &flag, MPI_STATUS_IGNORE));
    return flag != 0;
}

inline bool
Communicator::iprobe(int source, int tag, Status & status) const
{
    int flag;
    MPI_CHECK_SELF(MPI_Iprobe(source, tag, this->comm, &flag, status));
    return flag != 0;
}

// Barrier

inline void
Communicator::barrier() const
{
    MPI_CHECK_SELF(MPI_Barrier(this->comm));
}

// Broadcast

template <typename T>
void
Communicator::broadcast(T & value, int root) const
{
    broadcast(&value, 1, root);
}

template <typename T>
void
Communicator::broadcast(std::vector<T> & value, int root) const
{
    if (size() < 2)
        return;

    int tag = 0;
    if (rank() == root) {
        for (int i = 0; i < size(); ++i) {
            if (i != root)
                send(i, tag, value.data(), value.size());
        }
    }
    else
        recv(root, tag, value);
}

template <typename T>
void
Communicator::broadcast(T * values, int n, int root) const
{
    MPI_CHECK_SELF(MPI_Bcast(values, n, get_mpi_datatype<T>(), root, this->comm));
}

template <>
inline void
Communicator::broadcast(std::string & value, int root) const
{
    if (size() < 2)
        return;

    int tag = 0;
    if (rank() == root) {
        for (int i = 0; i < size(); ++i) {
            if (i != root)
                send(i, tag, value.data(), value.size());
        }
    }
    else {
        std::vector<char> str;
        recv(root, tag, str);
        value.assign(str.begin(), str.end());
    }
}

// Gather

template <typename T>
void
Communicator::gather(const T & in_value, T * out_values, int root) const
{
    assert(out_values || (rank() != root));
    gather(&in_value, 1, out_values, root);
}

template <typename T>
void
Communicator::gather(const T & in_value, std::vector<T> & out_values, int root) const
{
    if (rank() == root)
        out_values.resize(size());
    gather(in_value, out_values.data(), root);
}

template <typename T>
void
Communicator::gather(const T * in_values, int n, T * out_values, int root) const
{
    auto type = get_mpi_datatype<T>();
    MPI_CHECK_SELF(
        MPI_Gather(const_cast<T *>(in_values), n, type, out_values, n, type, root, this->comm));
}

template <typename T>
void
Communicator::gather(const T * in_values, int n, std::vector<T> & out_values, int root) const
{
    if (rank() == root)
        out_values.resize(size() * (std::size_t) n);
    gather(in_values, n, out_values.data(), root);
}

template <typename T>
void
Communicator::all_gather(const T * in_value, int n, T * out_values, int m) const
{
    MPI_CHECK_SELF(MPI_Allgather(in_value,
                                 n,
                                 get_mpi_datatype<T>(),
                                 out_values,
                                 m,
                                 get_mpi_datatype<T>(),
                                 this->comm));
}

template <typename T>
void
Communicator::all_gather(const T & in_value, std::vector<T> & out_values) const
{
    out_values.resize(this->size());
    all_gather(&in_value, 1, out_values.data(), 1);
}

template <typename T>
void
Communicator::all_gather(const std::vector<T> & in_values, std::vector<T> & out_values) const
{
    if (size() < 2)
        out_values = in_values;
    else {
        std::vector<int> n;
        int sz = in_values.size();
        all_gather(sz, n);
        std::vector<int> offsets(size());
        offsets[0] = 0;
        for (int i = 0; i < n.size() - 1; i++)
            offsets[i + 1] = offsets[i] + n[i];
        all_gather(in_values, out_values, n, offsets);
    }
}

template <typename T>
void
Communicator::all_gather(const std::vector<T> & in_values,
                         std::vector<T> & out_values,
                         const std::vector<int> & out_counts,
                         const std::vector<int> & out_offsets) const
{
    assert(out_counts.size() == size());
    assert(out_offsets.size() == size());
    int n_out_vals = 0;
    for (std::size_t i = 0; i < out_counts.size(); i++)
        n_out_vals += out_counts[i];
    out_values.resize(n_out_vals);
    MPI_CHECK_SELF(MPI_Allgatherv(in_values.data(),
                                  in_values.size(),
                                  get_mpi_datatype<T>(),
                                  out_values.data(),
                                  out_counts.data(),
                                  out_offsets.data(),
                                  get_mpi_datatype<T>(),
                                  this->comm));
}

// Scatter

template <typename T>
void
Communicator::scatter(const T * in_values, T & out_value, int root) const
{
    scatter(in_values, &out_value, 1, root);
}

template <typename T>
void
Communicator::scatter(const std::vector<T> & in_values, T & out_value, int root) const
{
    scatter(in_values.data(), &out_value, 1, root);
}

template <typename T>
void
Communicator::scatter(const T * in_values, T * out_values, int n, int root) const
{
    auto type = get_mpi_datatype<T>();
    MPI_CHECK_SELF(
        MPI_Scatter(const_cast<T *>(in_values), n, type, out_values, n, type, root, this->comm));
}

template <typename T>
void
Communicator::scatter(const std::vector<T> & in_values, T * out_values, int n, int root) const
{
    scatter(in_values.data(), out_values, n, root);
}

// Reduce

template <typename T, typename Op>
void
Communicator::reduce(const T * in_values, int n, T * out_values, Op, int root) const
{
    MPI_Op op = mpicpp_lite::op::Operation<Op, T>::op();
    MPI_CHECK_SELF(MPI_Reduce(const_cast<T *>(in_values),
                              out_values,
                              n,
                              mpicpp_lite::get_mpi_datatype<T>(),
                              op,
                              root,
                              this->comm));
}

template <typename T, typename Op>
void
Communicator::reduce(std::vector<T> const & in_values,
                     std::vector<T> & out_values,
                     Op op,
                     int root) const
{
    if (root == rank())
        out_values.resize(in_values.size());
    reduce(in_values.data(), in_values.size(), out_values.data(), op, root);
}

template <typename T, typename Op>
void
Communicator::reduce(const T & in_value, T & out_value, Op op, int root) const
{
    reduce(&in_value, 1, &out_value, op, root);
}

// All reduce

template <typename T, typename Op>
void
Communicator::all_reduce(const T * in_values, int n, T * out_values, Op) const
{
    MPI_Op op = mpicpp_lite::op::Operation<Op, T>::op();
    MPI_CHECK_SELF(MPI_Allreduce(const_cast<T *>(in_values),
                                 out_values,
                                 n,
                                 mpicpp_lite::get_mpi_datatype<T>(),
                                 op,
                                 this->comm));
}

template <typename T, typename Op>
void
Communicator::all_reduce(const T & in_value, T & out_value, Op op) const
{
    all_reduce(&in_value, 1, &out_value, op);
}

template <typename T, typename Op>
void
Communicator::all_reduce(T & in_value, Op op) const
{
    T out_value;
    all_reduce(&in_value, 1, &out_value, op);
    in_value = out_value;
}

template <typename T>
void
Communicator::all_to_all(const T * in_values, int n, T * out_values, int m) const
{
    MPI_CHECK_SELF(MPI_Alltoall(in_values,
                                n,
                                get_mpi_datatype<T>(),
                                out_values,
                                m,
                                get_mpi_datatype<T>(),
                                this->comm));
}

template <typename T>
void
Communicator::all_to_all(const std::vector<T> & in_values, std::vector<T> & out_values) const
{
    assert(in_values.size() == size());
    out_values.resize(size());
    all_to_all(in_values.data(), 1, out_values.data(), 1);
}

template <typename T>
void
Communicator::all_to_all(const std::vector<std::vector<T>> & in_values,
                         std::vector<T> & out_values) const
{
    std::vector<int> in_count(size(), 0);
    for (std::size_t i = 0; i < in_values.size(); i++)
        in_count[i] = in_values[i].size();
    std::vector<int> out_count(size(), 0);
    all_to_all(in_count, out_count);

    std::vector<int> in_offsets(size());
    in_offsets[0] = 0;
    for (std::size_t i = 0; i < in_count.size() - 1; i++)
        in_offsets[i + 1] = in_offsets[i] + in_count[i];
    std::size_t n_send_vals = 0;
    for (std::size_t i = 0; i < in_count.size(); i++)
        n_send_vals += in_count[i];

    std::vector<int> out_offsets(size());
    out_offsets[0] = 0;
    for (std::size_t i = 0; i < out_count.size() - 1; i++)
        out_offsets[i + 1] = out_offsets[i] + out_count[i];

    std::vector<T> in_buffer;
    in_buffer.reserve(n_send_vals);
    for (auto & vec : in_values)
        for (auto & v : vec)
            in_buffer.push_back(v);

    all_to_all(in_buffer, in_count, in_offsets, out_values, out_count, out_offsets);
}

template <typename T>
void
Communicator::all_to_all(const std::vector<T> & in_values,
                         const std::vector<int> & in_counts,
                         const std::vector<int> & in_offsets,
                         std::vector<T> & out_values,
                         const std::vector<int> & out_counts,
                         const std::vector<int> & out_offsets) const
{
    assert(in_counts.size() == size());
    assert(in_offsets.size() == size());
    assert(out_counts.size() == size());
    assert(out_offsets.size() == size());
    int n_receive_vals = 0;
    for (std::size_t i = 0; i < out_counts.size(); i++)
        n_receive_vals += out_counts[i];
    out_values.resize(n_receive_vals);
    MPI_CHECK_SELF(MPI_Alltoallv(in_values.data(),
                                 in_counts.data(),
                                 in_offsets.data(),
                                 get_mpi_datatype<T>(),
                                 out_values.data(),
                                 out_counts.data(),
                                 out_offsets.data(),
                                 get_mpi_datatype<T>(),
                                 this->comm));
}

//

inline void
Communicator::abort(int errcode) const
{
    MPI_Abort(this->comm, errcode);
}

inline Communicator::operator MPI_Comm() const
{
    return this->comm;
}

inline void
Communicator::set_error_handler()
{
#if (MPI_VERSION >= 2)
    MPI_Comm_set_errhandler(this->comm, MPI_ERRORS_RETURN);
#else
    MPI_Errhandler_set(this->comm, MPI_ERRORS_RETURN);
#endif
}

} // namespace mpicpp_lite
