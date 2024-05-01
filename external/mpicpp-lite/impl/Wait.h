// SPDX-FileCopyrightText: 2023 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpi.h"
#include <vector>
#include "Request.h"
#include "Status.h"
#include "Error.h"

namespace mpicpp_lite {

/// Wait for a single request to complete, ignoring status
///
/// @param request Request to wait for
inline void
wait(const Request & request)
{
    MPI_CHECK(MPI_Wait(const_cast<Request &>(request), MPI_STATUS_IGNORE));
}

/// Wait for a single request to complete with status
///
/// @param request Request to wait for
/// @param status Status
inline void
wait(const Request & request, Status & status)
{
    MPI_CHECK(MPI_Wait(const_cast<Request &>(request), status));
}

/// Wait for all requests to complete
///
/// @param requests Requests to wait for
inline void
wait_all(const std::vector<Request> & requests)
{
    int n = requests.size();
    std::vector<MPI_Request> reqs;
    reqs.resize(n);
    for (int i = 0; i < n; i++)
        reqs[i] = requests[i];
    MPI_CHECK(MPI_Waitall(n, reqs.data(), MPI_STATUSES_IGNORE));
}

/// Wait for any specified request to complete
///
/// @param requests Requests to wait for
/// @return Index of the request that completed
inline std::size_t
wait_any(const std::vector<Request> & requests)
{
    int n = requests.size();
    std::vector<MPI_Request> reqs;
    reqs.resize(n);
    for (int i = 0; i < n; i++)
        reqs[i] = requests[i];
    int idx;
    MPI_CHECK(MPI_Waitany(n, reqs.data(), &idx, MPI_STATUS_IGNORE));
    return idx;
}

} // namespace mpicpp_lite
