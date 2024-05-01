// SPDX-FileCopyrightText: 2023 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#ifdef MPICPP_LITE_WITH_FMT
    #include "fmt/printf.h"
#else
    #include <iostream>
#endif
#include "mpi.h"

namespace mpicpp_lite {

inline int
error_class(int code)
{
    int err_class;
    MPI_Error_class(code, &err_class);
    return err_class;
}

/// Return an error message for a given error code
///
/// @param code Error code
/// @return Error message
inline std::string
error_message(int code)
{
    char buffer[MPI_MAX_ERROR_STRING];
    int len;
    MPI_Error_string(code, buffer, &len);
    std::string msg(buffer);
    return msg;
}

namespace internal {

/// Terminate the run
inline void
terminate(MPI_Comm comm, int status = 1)
{
    MPI_Abort(comm, status);
}

inline void
check_mpi_error(MPI_Comm comm, int ierr, const char * file, int line)
{
    if (ierr) {
#ifdef MPICPP_LITE_WITH_FMT
        fmt::print(stderr,
                   "[ERROR] MPI error {} at {}:{}: {}\n",
                   ierr,
                   file,
                   line,
                   error_message(error_class(ierr)));
#else
        std::cerr << "[ERROR] MPI error " << ierr << " at " << file << ":" << line << ": "
                  << error_message(error_class(ierr)) << std::endl;
#endif
        terminate(comm, ierr);
    }
}

} // namespace internal

/// Check that MPI call was successful.
#define MPI_CHECK_SELF(ierr) \
    mpicpp_lite::internal::check_mpi_error(this->comm, ierr, __FILE__, __LINE__)

#define MPI_CHECK(ierr) \
    mpicpp_lite::internal::check_mpi_error(MPI_COMM_WORLD, ierr, __FILE__, __LINE__)

} // namespace mpicpp_lite
