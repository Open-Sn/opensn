/**\page DevManGlobalVars Global variables

\section devman2_sec0 Global variables available

All global entities for the ChiTech library are contained within the
static `chi` class. This class has the name `chi` and an instance of it
cannot be created. Within this class, all of its members are statically
declared. Several singleton objects are defined in `chi`, they are:
 - chi::mpi A handler for parallel related items.
 - chi::program_timer The primary program timer.
 - chi::console The link to the lua scripting engine.
 - chi::log A handler for parallel logging events and verbosity.

A number of stacks are also declared. They are basically arrays of
shared pointers (i.e., `std::shared_ptr`)

There are also a number of secondary global variables that assist developers
with coding. They are:
 - ```chi::run_time::input_file_name``` Holds the input file name if supplied.
 - ```chi::run_time::termination_posted``` A flag used during interactive mode.
 - ```chi::run_time::sim_option_interactive``` A flag indicating whether the code is
   run in interactive mode.
 - ```chi::run_time::allow_petsc_error_handler``` A flag indicating whether the allow
 the native PETSC error handler.



\subsection devman2_sec0_3 Connecting to the global data block

The stack items stored within `chi` can be accessed either by reference or as a
 `std::shared_ptr`. The two functions that facilitate this are as follows:

\code
#include "chi_runtime.h"
#include "mesh/SurfaceMesh/chi_surfacemesh.h" //Just an example

void SomeFunction()
{
  int handle = 2;
  auto& chi::GetStackItem<chi_mesh::SurfaceMeshPtr>(chi::surface_mesh_stack,
                                                    handle);
  //or
  auto chi::GetStackItemPtr<chi_mesh::SurfaceMeshPtr>(chi::surface_mesh_stack,
                                                      handle);
}
\endcode

There are multiple stacks, currently (which might not be still accurate):
\code
chi::meshhandler_stack;
chi::surface_mesh_stack;
chi::logicvolume_stack;
chi::field_func_interpolation_stack;
chi::unpartitionedmesh_stack;
chi::solver_stack;
chi::material_stack;
chi::trnsprt_xs_stack;
chi::fieldfunc_stack;
chi::quadrature_stack;
chi::angular_quadrature_stack;
\endcode

\subsection devman2_sec0_4 Connecting to MPI

General MPI information like the current location id and the total amount
 of parallel processes is contained in CHI_MPI:
 - objects::MPI_Info::location_id
 - objects::MPI_Info::process_count

Additionally, by including the headers for chi_mpi, developers have access to
 all the classic mpi headers.

\code
#include "chi_runtime.h"
#include "chi_mpi.h"

if (chi::mpi.location_id == 1)
 printf("This is process 1, Dab!");
\endcode

\subsection devman2_sec0_5 Connecting to the parallel logging utility

Printing information in a parallel environment can be a very involved
process. One can't simply use `std::cout <<` on every process otherwise
the output to the log will be chaotic. For this reason we employ a common
 logging utility which returns an output string-stream using the function
 call Log::Log.

Connecting to chi::log is done as follows
\code
#include "chi_runtime.h"
#include "chi_log.h"

void Function()
{
    chi::log.Log() << "Hello from location 0";
    chi::log.LogAll() << "Hello from all locations";
}
\endcode

The logger has calls of differing verbosity:
 - `objects::Log::Log0()`,                      Used only for location 0
 - `objects::Log::Log0Warning()`,               Warning only for location 0
 - `objects::Log::Log0Error()`,                 Error only for location 0
 - `objects::Log::Log0Verbose0()`,             Default verbosity level
 - `objects::Log::Log0Verbose1()`,             Used only if verbosity level equals 1
 - `objects::Log::Log0Verbose2()`,             Used only if verbosity level equals 2
 - `objects::Log::LogAll()`,                    Verbose level 0 all locations
 - `objects::Log::LogAllWarning()`,             Warning for any location
 - `objects::Log::LogAllError()`,               Error for any location
 - `objects::Log::LogAllVerbose0()`,     Default verbosity level
 - `objects::Log::LogAllVerbose1()`,     Used only if verbosity level equals 1
 - `objects::Log::LogAllVerbose2()`,     Used only if verbosity level equals 2



\todo Commandline parameter documentation


*/
