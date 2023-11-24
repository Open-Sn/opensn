#pragma once

#include "mpi.h"

namespace mpicpp_lite {

/// Wrapper around MPI_Request
class Request {
public:
    /// Create empty request
    Request();

    /// Create request from an `MPI_Request`
    ///
    /// @param r `MPI_Request` used to initiliaze this object
    Request(const MPI_Request & r);

    /// Type cast operator so we can pass this class directly into MPI calls
    operator MPI_Request *() { return &this->request; }

    /// Type cast operator so we can pass this class directly into MPI calls
    operator const MPI_Request &() const { return this->request; }

private:
    MPI_Request request;
};

inline Request::Request() {}

inline Request::Request(const MPI_Request & r) : request(r) {}

} // namespace mpi
