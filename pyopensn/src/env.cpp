// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include "mpi4py/mpi4py.h"
#include "petscsys.h"

#include "framework/runtime.h"
#include "framework/utils/timer.h"

namespace opensn {

// Constructor
PyEnv::PyEnv(void) {
    // check if environment is already initialized
    if (PyEnv::p_default_env != nullptr) {
        return;
    }
    // import mpi4py and get the communicator
    ::import_mpi4py();
    py::object mpi_module = py::module::import("mpi4py.MPI");
    py::object comm_world = mpi_module.attr("COMM_WORLD");
    mpi_comm = mpi::Communicator(*(::PyMPIComm_Get(comm_world.ptr())));
    if (mpi_comm.rank() == 0) {
        std::cout << program << " version " << GetVersionStr() << "\n";
        std::cout << Timer::GetLocalDateTimeString() << " Running " << program << " with " << mpi_comm.size()
                  << " processes.\n\n";
    }
    // initialize PETSc and Opensn
    ::PetscOptionsSetValue(NULL, "-options_left", "0");
    ::PetscOptionsInsertString(nullptr, "-no_signal_handler");
    ::PetscInitialize(nullptr, nullptr, nullptr, nullptr);
    Initialize();
}

// Destructor
PyEnv::~PyEnv(void) {
    // finalize the run
    Finalize();
    ::PetscFinalize();
    // print execution time
    if (mpi_comm.rank() == 0) {
        std::cout << "\n";
        std::cout << "Elapsed execution time: " << program_timer.GetTimeString() << "\n";
        std::cout << Timer::GetLocalDateTimeString() << " " << program << " finished execution.\n";
    }
    // flush caliper
    cali_mgr.flush();
}

// Default environement
PyEnv * PyEnv::p_default_env = nullptr;

}  // namespace opensn
