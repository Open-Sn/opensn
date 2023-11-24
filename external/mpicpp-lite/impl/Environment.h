#pragma once

#include "mpi.h"
#include "Error.h"

namespace mpicpp_lite {

class Environment {
public:
#if (MPI_VERSION >= 2)
    Environment();
#endif
    Environment(int& argc, char** &argv);
    ~Environment();

    /// Indicates whether MPI was initialized
    ///
    /// @return `true` if MPI is initialized, `false` otherwise
    static bool is_initialized();

    /// Indicates whether MPI was finalized
    ///
    /// @return `true` if MPI is finalized, `false` otherwise
    static bool is_finalized();

private:
    /// Indicates if the environment is initialized
    bool initialized;
};

inline
Environment::Environment() : initialized(false)
{
    if (!is_initialized()) {
        MPI_CHECK(MPI_Init(nullptr, nullptr));
        this->initialized = true;
    }
}

inline
Environment::Environment(int& argc, char** &argv) : initialized(false)
{
    if (!is_initialized()) {
        MPI_CHECK(MPI_Init(&argc, &argv));
        this->initialized = true;
    }
}

inline
Environment::~Environment()
{
    if (this->initialized) {
        if (!is_finalized())
            MPI_CHECK(MPI_Finalize());
    }
}

inline bool
Environment::is_initialized()
{
    int flag;
    MPI_CHECK(MPI_Initialized(&flag));
    return flag != 0;
}

inline bool
Environment::is_finalized()
{
    int flag;
    MPI_CHECK(MPI_Finalized(&flag));
    return flag != 0;
}

} // namespace mpicpp_lite
