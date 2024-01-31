#pragma once

#include "mpi.h"
#include <vector>
#include "Request.h"
#include "Status.h"
#include "Error.h"

namespace mpicpp_lite {

/// Test for the completion of a request
///
/// @param request Request to test
/// @return `true` if operation completed, `false` otherwise
inline bool
test(const Request & request)
{
    int flag;
    MPI_CHECK(MPI_Test(const_cast<Request &>(request), &flag, MPI_STATUS_IGNORE));
    return flag != 0;
}

/// Test for the completion of a request with status
///
/// @param request Request to wait for
/// @param status Status
/// @return `true` if operation completed, `false` otherwise
inline bool
test(const Request & request, Status & status)
{
    int flag;
    MPI_CHECK(MPI_Test(const_cast<Request &>(request), &flag, status));
    return flag != 0;
}

/// Test for the completion of all previously initiated requests
///
/// @param requests Requests to test
/// @return `true` only if all requests have completed, `false` otherwise
inline bool
test_all(const std::vector<Request> & requests)
{
    int n = requests.size();
    std::vector<MPI_Request> reqs;
    reqs.resize(n);
    for (int i = 0; i < n; i++)
        reqs[i] = requests[i];
    int flag;
    MPI_CHECK(MPI_Testall(n, reqs.data(), &flag, MPI_STATUSES_IGNORE));
    return flag != 0;
}

/// Test for completion of any previously initiated requests
///
/// @param requests Requests to test
/// @param index Index of operation that completed or `UNDEFINED` if none completed
/// @return `true` if one of the operations is complete
inline bool
test_any(const std::vector<Request> & requests, std::size_t & index)
{
    int n = requests.size();
    std::vector<MPI_Request> reqs;
    reqs.resize(n);
    for (int i = 0; i < n; i++)
        reqs[i] = requests[i];
    int idx;
    int flag;
    MPI_CHECK(MPI_Testany(n, reqs.data(), &idx, &flag, MPI_STATUS_IGNORE));
    return flag != 0;
}

} // namespace mpicpp_lite
