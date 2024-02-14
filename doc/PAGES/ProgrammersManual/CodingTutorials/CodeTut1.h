/** \page CodeTut1 Coding Tutorial 1 - Poisson's problem with Finite Volume

## Table of contents
- \ref CodeTut1Sec1
 - \ref CodeTut1Sec1_1
 - \ref CodeTut1Sec1_2
- \ref CodeTut1Sec2
- \ref CodeTut1Sec3
- \ref CodeTut1Sec4
- \ref CodeTut1Sec5
- \ref CodeTut1Sec6
- \ref CodeTut1Sec7
- \ref CodeTut1Sec8
- \ref CodeTut1Sec9

\section CodeTut1Sec1 1 Poisson's equation
<a href="https://en.wikipedia.org/wiki/Poisson%27s_equation">The Poisson's equation</a> states the
following, for \f$ \phi, q \in \mathcal{R} \f$, \f{eqnarray*}{
-\boldsymbol{\nabla} \cdot \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) = q(\mathbf{x}), \quad
\quad \quad (\mathbf{x})\in\mathcal{D} \\ \phi(\mathbf{x}) = 0, \quad \quad \quad \mathbf{x}\in
\partial \mathcal{D} \f} where \f$ \boldsymbol{\nabla} \cdot \bigr( \bigr) \f$ denotes the
divergence-operator and \f$ \boldsymbol{\nabla} \f$ denotes the gradient-operator, i.e. \f$
\boldsymbol{\nabla} \phi \f$ denotes the gradient of \f$ \phi \f$. The boundary conditions here
state that \f$ \phi=0 \f$ on the boundary.




\subsection CodeTut1Sec1_1 1.1 Our specific problem
For our specific problem we will choose \f$ q(\mathbf{x})=1 \f$ and \f$ \mathcal{D} \f$ a cartesian
domain, either 1D, 2D or 3D, with each dimension always between \f$ -1,+1 \f$. We can generate the
mesh for this problem using an input file \code
--############################################### Setup mesh
MeshHandlerCreate()

mesh={}
N=20
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

--MeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
MeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--MeshCreateUnpartitioned1DOrthoMesh(mesh)
VolumeMesherExecute();

--############################################### Set Material IDs
mesh.SetUniformMaterialID(0)
\endcode
This code can be used to generate any of the following meshes,
\image html CodingTutorials/OrthoMesh_1D_2D_3D.png "[From left to right] 1D, 2D, 3D orthogonal mesh"
width=1200px




\subsection CodeTut1Sec1_2 1.2 The Finite Volume spatial discretization
To apply the Finite Volume spatial discretization we integrate Poisson's equation
over the volume of the cell
\f{eqnarray*}{
-\int_V \boldsymbol{\nabla} \cdot \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) dV
=
\int_V q(\mathbf{x}) dV
\f}
Next we apply the approximation that both \f$ \phi \f$ and \f$ q \f$ are constant
over a cell. Therefore, for cell \f$ P \f$, we have the source term
\f{eqnarray*}{
\int_V q(\mathbf{x}) dV = q_P V_P.
\f}
For the left-hand side we first apply Gauss Divergence Theorem,
\f{eqnarray*}{
\int_V \boldsymbol{\nabla} \cdot \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) dV =
\int_S \mathbf{n} \cdot \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) dA
\f}
and for a discrete cell this surface integral can be discretized as the sum over all
the faces
\f{eqnarray*}{
\int_S \mathbf{n} \cdot \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr) dA =
\sum_f A_f \mathbf{n}_f \cdot \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr)_f.
\f}
We now still have to resolve \f$ \bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr)_f \f$.
To comprehend the approximation we are about to make, consider the figure below.
We observe a cell \f$ P \f$, the present cell, and its neighbor, cell \f$ N \f$, at face \f$ f \f$.
The cell centroids are \f$ \mathbf{x}_P \f$  and \f$ \mathbf{x}_N \f$ respectively.

\image html CodingTutorials/NeighborXPN.png width=250px

We next form \f$ \mathbf{x}_{PN} \f$, as in the diagram, after which we approximate
the gradient of \f$ \phi \f$ at face \f$ f \f$ as
\f{eqnarray*}{
\bigr(\boldsymbol{\nabla} \phi(\mathbf{x})\bigr)_f =
\biggr(
\frac{\phi_N - \phi_P}{||\mathbf{x}_{PN}||}
\biggr)
\ \frac{\mathbf{x}_{PN}}{||\mathbf{x}_{PN}||}.
\f}
Note that this is essentially a scalar component
 \f$ \frac{df}{ds} = \frac{\phi_N - \phi_P}{||\mathbf{x}_{PN}||}\f$ after which we give it a
direction by normalizing \f$ \mathbf{x}_{PN} \f$ as \f$ \frac{\mathbf{x}_{PN}}{||\mathbf{x}_{PN}||}
\f$.

Our discretized equations now become

\f{align*}{
\sum_f
c_f
\phi_P
-
\sum_f
c_f
\phi_N
&=
q_P V_P
\\
c_f &= \biggr[
\mathbf{A}_f \boldsymbol{\cdot}
\frac{\mathbf{x}_{PN}}{||\mathbf{x}_{PN}||^2}
\biggr]_f
\f}
where \f$ \mathbf{A}_f \f$ is the area-vector of the face, i.e. \f$ A_f \mathbf{n}_f \f$ and
 \f$ \phi_P=0 \f$ when face \f$ f \f$ is on a boundary (enforcing the
 Dirichlet boundary condition).





\section CodeTut1Sec2 2 Setting up the problem
In order to implement the code, we must first follow the steps in \ref DevManUsingLib.
Before proceeding, make sure you are fully acquinted with this step.

For this tutorial we will be deviation from the library tutorial by first
changing the cmake target and project name. Therefore, we pick a project folder
and create the `CMakeLists.txt` file but we change the following lines.
\code
set(TARGET code_tut1)
project(CodeTut1 CXX)
\endcode

So now the entire `CMakeLists.txt` should look like this
\code
cmake_minimum_required(VERSION 3.12)

set(TARGET code_tut1)
project(CodeTut1 CXX)

#------------------------------------------------ DEPENDENCIES
if (NOT DEFINED CHI_TECH_DIR)
    if (NOT (DEFINED ENV{CHI_TECH_DIR}))
        message(FATAL_ERROR "***** CHI_TECH_DIR is not set *****")
    else()
        set(CHI_TECH_DIR "$ENV{CHI_TECH_DIR}")
    endif()
endif()
message(STATUS "CHI_TECH_DIR set to ${CHI_TECH_DIR}")

include("${CHI_TECH_DIR}/resources/Macros/Downstream.cmake")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

file (GLOB_RECURSE SOURCES "*.cc")
add_executable(${TARGET} "${SOURCES}")
target_link_libraries(${TARGET} ${CHI_LIBS})

file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build \n\n"
        "clean:\n\t$(MAKE) -C chi_build clean\n")
\endcode

Next we create the source file which we will call `code_tut1.cc`, and for now
we will keep its contents the same as the library tutorial. Therefore `code_tut1.cc`
should look like this
\code
#include "chi_runtime.h"
#include "chi_log.h"

int main(int argc, char* argv[])
{
  chi::Initialize(argc,argv);
  chi::RunBatch(argc, argv);

  chi::log.Log() << "Hello World!";

  //We will add code here

  chi::Finalize();
  return 0;
}
\endcode

To compile this program we first have to set the `CHI_TECH_DIR` environment
variable. This might be different for everyone but should generally look like this
\code
export CHI_TECH_DIR=/Users/drjanisawesome/codes/ChiTech/chi-tech/
\endcode
Note: the directory we want to specify here must contain `bin/` (in otherwords)
it shouldn't be `bin/` itself).

Next we create a directory for running `cmake`. Inside your project folder, create
a folder called `chi_build` and `cd` into it.
\code
mkdir chi_build
cd chi_build
\endcode

Then run `cmake` pointing to the folder that contains the `CMakeLists.txt` file,
which in this case is the parent folder of the one we are in.
\code
cmake ../
\endcode

After this you should be able to make the project.
\code
make
\endcode
after which you should see the following
\verbatim
[ 50%] Building CXX object CMakeFiles/code_tut1.dir/code_tut1.cc.o
[100%] Linking CXX executable ../bin/code_tut1
[100%] Built target code_tut1
\endverbatim

As a test you can try to execute the project with no arguments which should look
like this
\verbatim
../bin/code_tut1
[0]  2022-11-14 13:02:43 Running ChiTech in batch-mode with 1 processes.
[0]  ChiTech version 1.0.2
[0]  ChiTech number of arguments supplied: 0
[0]
[0]  Usage: exe inputfile [options values]
[0]
[0]       -v                         Level of verbosity. Default 0. Can be either 0, 1 or 2.
[0]       a=b                        Executes argument as a lua string. i.e. x=2 or y=[["string"]]
[0]       -allow_petsc_error_handler Allow petsc error handler.
[0]
[0]
[0]
[0]  Final program time 00:00:00
[0]  2022-11-14 13:02:43 ChiTech finished execution of
[0]  Hello World!
\endverbatim

\section CodeTut1Sec3 3 Getting the grid
First we remove the lines
\code
chi::log.Log() << "Hello World!";

//We will add code here
\endcode

Next we include the headers for access to the `chi_mesh::MeshHandler` and the
grid, aka `chi_mesh::MeshContinuum`
\code
#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
\endcode

Next we add the following lines of code:
\code
chi::log.Log() << "Coding Tutorial 1";

// Get grid
auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
const auto& grid = *grid_ptr;

chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();
\endcode
This is the default way of obtaining a grid loaded into the ChiTech state at any
moment, normally via an input file.

When you compile and run this code it ought to fail by throwing an exception
that looks like this
\verbatim
terminate called after throwing an instance of 'std::logic_error'
  what():  chi_mesh::GetCurrentHandler: No handlers on stack
Abort trap: 6
\endverbatim
This means the code was not loaded with a mesh-handler and therefore no grid is
defined. To fix this we first create a `mesh.lua` input file within the `chi_build`
directory.
\code
>mesh.lua
\endcode
Use an editor to add the following lines to `mesh.lua`
\code
--############################################### Setup mesh
MeshHandlerCreate()

mesh={}
N=20
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end

--MeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
MeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--MeshCreateUnpartitioned1DOrthoMesh(mesh)
VolumeMesherExecute();

--############################################### Set Material IDs
mesh.SetUniformMaterialID(0)
\endcode
This is known as ChiTech lua-console input file. It contains lua code.

To run the program with this input file we simply add it to the command line
\code
../bin/code_tut1 mesh.lua
\endcode

If everything went well then you should see the following
\verbatim
[0]  Parsing argument 1 mesh.lua
[0]  2022-11-14 13:14:52 Running ChiTech in batch-mode with 1 processes.
[0]  ChiTech version 1.0.2
[0]  ChiTech number of arguments supplied: 1
[0]  Computing cell-centroids.
[0]  Done computing cell-centroids.
[0]  Checking cell-center-to-face orientations
[0]  Done checking cell-center-to-face orientations
[0]  00:00:00 Establishing cell connectivity.
[0]  00:00:00 Vertex cell subscriptions complete.
[0]  00:00:00 Surpassing cell 40 of 400 (10%)
[0]  00:00:00 Surpassing cell 80 of 400 (20%)
\\ STUFF
\\ STUFF
\\ STUFF
[0]
[0]  Final program time 00:00:00
[0]  2022-11-14 13:14:52 ChiTech finished execution of mesh.lua
[0]  Coding Tutorial 1
[0]  Global num cells: 400
\endverbatim

\section CodeTut1Sec4 4 Initializing the Finite Volume Spatial Discretization
All spatial discretizations in ChiTech are based on a grid. Right now we want
`chi_math::SpatialDiscretization_FV` which is the Finite Volume Spatial discretization.
It will place only a single node at the centroid of each cell.

To access this discretization we must first include the header for
`chi_math::SpatialDiscretization_FV`,
\code
#include "math/SpatialDiscretization/FiniteVolume/fv.h"
\endcode

Next we add the following lines of code
\code
// Make SDM
typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
SDMPtr sdm_ptr = chi_math::SpatialDiscretization_FV::New(grid_ptr);
const auto& sdm = *sdm_ptr;

const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

chi::log.Log() << "Num local DOFs: " << num_local_dofs;
chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;
\endcode

What we did here is to first define the shared-pointer of type
`chi_math::SpatialDiscretization` to be just `SMDPtr` which greatly reduced
the effort on the next line, where we actually instantiate the spatial discretization.
Notice the use of `chi_math::SpatialDiscretization_FV::New`. All of the ChiTech
spatial discretizations can only be created as shared pointers, the need for this
will only become obvious when you write your own solvers.

The creation of a spatial discretization involves a lot of steps. The first thing
the spatial discretization (SD) will do is to create cell-mappings. These contain
all the necessary information for each cell to help build codes, e.g., the number
of nodes on a cell, its volume, face areas, etc. The second thing the SD does is
to order the nodes in parallel in such a way that we can easily map the nodes either
globally or locally.

The basic SD operates on the notion of `nodes`. A node is a place where a specific
component of an `unknown` is located. Since it knows the underlying node structure
it can provide index mappings for different unknown-structures. For example, if we
stack velocity-x, velocity-y and velocity-z onto a node the SD only needs to know
the structure of the unknowns, traditionally called Degrees-of-Freedom or DOFs, in
order to provide index-mapping. To keep things simple, all SDs have a unitary
`chi_math::UnknownManager` as a member called `UNITARY_UNKNOWN_MANAGER`. This
unknown manager allows us to map indices assuming only one DOF per node. We obtained
a reference to this `chi_math::UnknownManager` object with the line
\code
const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
\endcode

With this in-hand we can now query the SD to provide us the number of local- and
global DOFs. The local DOFs are those that are stored on the local processor whereas
the global DOFs encompass all the DOFs across all processors. The way we ran the
code above is only in serial. Therefore when we do this again, with this new code
added, we will see the following new output:
\verbatim
[0]  Num local DOFs: 400
[0]  Num globl DOFs: 400
\endverbatim
Since the code runs in serial mode the number of local- and global DOFS are reported
to be the same.

To run this code in parallel we change our input, for execution with 3 processors,
as
\code
mpiexec -np 3 ../bin/code_tut1 mesh.lua
\endcode
which now produces a different output:
\verbatim
[0]  Num local DOFs: 131
[0]  Num globl DOFs: 400
\endverbatim
Notice the global number of DOFs remains the same. Only the number of local nodes
has changed.

\section CodeTut1Sec5 5 Creating and Initializing PETSc matrices and vectors
ChiTech has several `macro`-type functions for handling PETSc objects. All of these
are accessed via the `chi_math::PETScUtils` namespace for which we need to include
the header
\code
#include "math/PETScUtils/petsc_utils.h"
\endcode

Next we add the following code:
\code
// Initializes Mats and Vecs
const auto n = static_cast<int64_t>(num_local_dofs);
const auto N = static_cast<int64_t>(num_globl_dofs);
Mat A;
Vec x,b;

A = chi_math::PETScUtils::CreateSquareMatrix(n,N);
x = chi_math::PETScUtils::CreateVector(n,N);
b = chi_math::PETScUtils::CreateVector(n,N);

std::vector<int64_t> nodal_nnz_in_diag;
std::vector<int64_t> nodal_nnz_off_diag;
sdm.BuildSparsityPattern(nodal_nnz_in_diag,nodal_nnz_off_diag,OneDofPerNode);

chi_math::PETScUtils::InitMatrixSparsity(A,
                                         nodal_nnz_in_diag,
                                         nodal_nnz_off_diag);
\endcode
Notice that we made `int64_t` equivalents of `num_local_dofs` and its global
counterpart. This is because PETSc operates on `int64_t`. The actual creation of
`A`, `x`, and `b` is self explanatory noting that each needs a local-size as well
as a global size.

The next two pieces of code after the creation of the PETSc entities allow us to
accurately and efficiently set the matrix sparsity pattern, in-turn allowing PETSc to
assemble the sparse-matrix efficiently. The sparsity pattern is defined per row
and is split into entries locally stored ("in" the diagonal block) and entries
stored elsewhere ("off" the diagonal block). We can get this information from our
ChiTech SDs by using the `BuildSparsityPattern` routine. This routine populates
`nodal_nnz_in_diag` and `nodal_nnz_off_diag`, according to an unknown-structure
(which in our case is `OneDofPerNode`).

Finally we use our sparsity information to set the PETSc matrix's sparsity
pattern using `chi_math::PETScUtils::InitMatrixSparsity`.






\section CodeTut1Sec6 6 Assembling the matrix and right-hand-side
The code to assemble the matrix is as follows.
\code
// Assemble the system
chi::log.Log() << "Assembling system: ";
for (const auto& cell : grid.local_cells)
{
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const int64_t imap = sdm.MapDOF(cell,0);

  const auto& xp = cell.centroid;
  const double V = cell_mapping.CellVolume();

  size_t f=0;
  for (const auto& face : cell.faces)
  {
    const auto Af = face.normal * cell_mapping.FaceArea(f);

    if (face.has_neighbor)
    {
      const auto& adj_cell = grid.cells[face.neighbor_id];
      const int64_t jnmap = sdm.MapDOF(adj_cell,0);

      const auto& xn = adj_cell.centroid;

      const auto xpn = xn - xp;

      const auto cf = Af.Dot(xpn) / xpn.NormSquare();

      MatSetValue(A, imap, imap ,  cf, ADD_VALUES);
      MatSetValue(A, imap, jnmap, -cf, ADD_VALUES);
    }
    else
    {
      const auto& xn = xp + 2.0*(face.centroid - xp);
      const auto xpn = xn - xp;

      const auto cf = Af.Dot(xpn) / xpn.NormSquare();

      MatSetValue(A, imap, imap , cf, ADD_VALUES);
    }
    ++f;
  }//for face

  VecSetValue(b, imap, 1.0*V, ADD_VALUES);
}//for cell i

chi::log.Log() << "Global assembly";

MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(b);
VecAssemblyEnd(b);

chi::log.Log() << "Done global assembly";
\endcode
The first item to note is the loop over local cells. ChiTech can only loop over local
cells, for more information on this see \ref devman_meshes_sec_1.

Once we have a reference to a local-cell we can now obtain the SD's mapping of that cell
using
\code
const auto& cell_mapping = sdm.GetCellMapping(cell);
\endcode
We can also immediately map the cell's \f$ \phi_P \f$ global location in the
parallel matrix/vector by making a call to
`chi_math::spatial_discretization::MapDOF`.
\code
const int64_t imap = sdm.MapDOF(cell,0);
\endcode
The syntax of this function is #1 A cell reference, and #2 the node number which,
for Finite Volume, will always be zero.

Finally, before looping over the faces, we grab the cell centroid and volume.

The loop over faces is split into two cases, one for internal faces and one for
boundary faces. Both cases need the face area-vector there we construct that first
\code
for (const auto& face : cell.faces)
{
  const auto Af = face.normal * cell_mapping.FaceArea(f);
  //etc.
\endcode

For internal faces, whether that face is internal, i.e., has a neighboring cell,
is indicated by the member `has_neighbor`. Additionally the neighbor's global-id
will be stored in the member `neighbor_id`. If the face does not have a neighbor
then `neighbor_id` doubles as the boundary-id if boundary-ids are used. All cells
that at minimum share a vertex with local cells are classified as `Ghost`-cells
in ChiTech. And since a neighbor can be either a local- or ghost-cell we opt here
to get a reference to the neighbor cell using its global-id for which we use
`grid.cells`, instead of `grid.local_cells` which only uses local-ids.
\code
const auto& adj_cell = grid.cells[face.neighbor_id];
\endcode
Once we have a reference to the adjacent cell we can now map the address of its
unknown \f$ \phi_N \f$ as
\code
const int64_t jnmap = sdm.MapDOF(adj_cell,0);
\endcode
We then compute
\f{align*}{
c_f = \biggr[
\mathbf{A}_f \boldsymbol{\cdot}
\frac{\mathbf{x}_{PN}}{||\mathbf{x}_{PN}||^2}
\biggr]_f
\f}
with the code
\code
const auto& xn = adj_cell.centroid;

const auto xpn = xn - xp;

const auto cf = Af.Dot(xpn) / xpn.NormSquare();
\endcode
after that we set the matrix entries associated with \f$ \phi_P \f$ and
\f$ \phi_N \f$ using the PETSc commands
\code
MatSetValue(A, imap, imap ,  cf, ADD_VALUES);
MatSetValue(A, imap, jnmap, -cf, ADD_VALUES);
\endcode




For boundary faces, we essentially do the same, however, we do not have an adjacent
cell or \f$ \mathbf{x}_N \f$. Therefore we extrapolate a virtual boundary cell-centroid
using the face-centroid, \f$ \mathbf{x}_f \f$, as
\f[
\mathbf{x}_N = \mathbf{x}_P + 2 ( \mathbf{x}_f - \mathbf{x}_P )
\f]
and then
\code
const auto& xn = xp + 2.0*(face.centroid - xp);
const auto xpn = xn - xp;

const auto cf = Af.Dot(xpn) / xpn.NormSquare();
\endcode

Finally, since there is no neighbor cell entry, we only have a diagonal entry
\code
MatSetValue(A, imap, imap , cf, ADD_VALUES);
\endcode

After the face loops we then set the right-hand side using
\code
VecSetValue(b, imap, 1.0*V, ADD_VALUES);
\endcode
where the `1.0` represents \f$ q(\mathbf{x})=1 \f$.

We finish parallel assembly with the following calls
\code
MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
VecAssemblyBegin(b);
VecAssemblyEnd(b);
\endcode




\section CodeTut1Sec7 7 Solving the system with a Krylov Subspace solver
Most of the following code is self explanatory.
\code
// Create Krylov Solver
chi::log.Log() << "Solving: ";
auto petsc_solver =
  chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
    A,               //Matrix
    "FVDiffSolver",  //Solver name
    KSPCG,           //Solver type
    PCGAMG,          //Preconditioner type
    1.0e-6,          //Relative residual tolerance
    1000);            //Max iterations

// Solve
KSPSolve(petsc_solver.ksp,b,x);

chi::log.Log() << "Done solving";
\endcode
The solver type here was specified as `KSPCG` which represents the conjugate-gradient
method. The preconditioner `PCGAMG` represents Geometric and/or Algebraic MultiGrid.
The combination of this solver and preconditioner is known to be thee most efficient
ways to solve these type of systems.

After this solve-step we want to get an STL representation of the PETSc vector so
that we can think about visualizing the solution. We can get a local copy of the
PETSc vector using the code
\code
// Extract PETSc vector
std::vector<double> field(num_local_dofs, 0.0);
sdm.LocalizePETScVector(x,field,OneDofPerNode);

// Clean up
KSPDestroy(&petsc_solver.ksp);

VecDestroy(&x);
VecDestroy(&b);
MatDestroy(&A);

chi::log.Log() << "Done cleanup";
\endcode
Here the SD took care of appropriately copying from the PETSc vector.

After this is completed we have no need for any of the PETSc entities and therefore
we destroy them.






\section CodeTut1Sec8 8 Visualizing the solution
ChiTech uses the notion of `FieldFunctions` and the visualization toolkit, `VTK`,
to visual solutions. In order to gain access to `chi_physics::FieldFunction` we
need to include the header
\code
#include "physics/FieldFunction/fieldfunction2.h"
\endcode

Next we create the field function using the code
\code
// Create Field Function
auto ff = std::make_shared<chi_physics::FieldFunction>(
  "Phi",
  sdm_ptr,
  chi_math::Unknown(chi_math::UnknownType::SCALAR)
);
\endcode
Note here that a field function simply needs a name, a SD, and an unknown structure.

After this is created we can update the field function with a field vector
\code
ff->UpdateFieldVector(field);
\endcode

And finally we can export the solution to VTK with
\code
ff->ExportToVTK("CodeTut1_FV");
\endcode
which takes, as an argument, the base-name of the VTK output files. VTK will
utilize the parallel nature of the simulation to generate a set of files
\verbatim
CodeTut1_FV.pvtu
CodeTut1_FV_0.vtu
CodeTut1_FV_1.vtu
CodeTut1_FV_2.vtu
\endverbatim
The structure here is that each processor writes a `.vtu` file. The processor-id
is appended to the basename before the `.vtu` extension is added. These files contain
the actual data. An additional proxy-file, used to link all the data together, is written
as the base-name with `.pvtu` appended to it. This `.pvtu` file can be opened in popular
visualization tools such as `Visit` and `Paraview` to visualize the solution.

For this tutorial we used Paraview.

\image html CodingTutorials/Solution_1D_2D_3D.png "[From left to right] 1D, 2D, 3D solutions with a
coarse mesh" width=1200px

And on a finer mesh, with the 3D case being 1 million cells:
\image html CodingTutorials/Solution_1D_2D_3D_fine.png "[From left to right] 1D, 2D, 3D solutions
with a fine mesh" width=1200px




\section CodeTut1Sec9 The complete program

\code
#include "chi_runtime.h"
#include "chi_log.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteVolume/fv.h"
#include "math/PETScUtils/petsc_utils.h"

#include "physics/FieldFunction/fieldfunction2.h"

int main(int argc, char* argv[])
{
  chi::Initialize(argc,argv);
  chi::RunBatch(argc, argv);

  chi::log.Log() << "Coding Tutorial 1";

  // Get grid
  auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr;

  chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make SDM
  typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
  SDMPtr sdm_ptr = chi_math::SpatialDiscretization_FV::New(grid_ptr);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

  chi::log.Log() << "Num local DOFs: " << num_local_dofs;
  chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs);
  const auto N = static_cast<int64_t>(num_globl_dofs);
  Mat A;
  Vec x,b;

  A = chi_math::PETScUtils::CreateSquareMatrix(n,N);
  x = chi_math::PETScUtils::CreateVector(n,N);
  b = chi_math::PETScUtils::CreateVector(n,N);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag,nodal_nnz_off_diag,OneDofPerNode);

  chi_math::PETScUtils::InitMatrixSparsity(A,
                                           nodal_nnz_in_diag,
                                           nodal_nnz_off_diag);

  // Assemble the system
  chi::log.Log() << "Assembling system: ";
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const int64_t imap = sdm.MapDOF(cell,0);

    const auto& xp = cell.centroid;
    const double V = cell_mapping.CellVolume();

    size_t f=0;
    for (const auto& face : cell.faces)
    {
      const auto Af = face.normal * cell_mapping.FaceArea(f);

      if (face.has_neighbor)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id];
        const int64_t jnmap = sdm.MapDOF(adj_cell,0);

        const auto& xn = adj_cell.centroid;

        const auto xpn = xn - xp;

        const auto cf = Af.Dot(xpn) / xpn.NormSquare();

        MatSetValue(A, imap, imap ,  cf, ADD_VALUES);
        MatSetValue(A, imap, jnmap, -cf, ADD_VALUES);
      }
      else
      {
        const auto& xn = xp + 2.0*(face.centroid - xp);
        const auto xpn = xn - xp;

        const auto cf = Af.Dot(xpn) / xpn.NormSquare();

        MatSetValue(A, imap, imap , cf, ADD_VALUES);
      }
      ++f;
    }//for face

    VecSetValue(b, imap, 1.0*V, ADD_VALUES);
  }//for cell i

  chi::log.Log() << "Global assembly";

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  chi::log.Log() << "Done global assembly";

  // Create Krylov Solver
  chi::log.Log() << "Solving: ";
  auto petsc_solver =
    chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
      A,               //Matrix
      "FVDiffSolver",  //Solver name
      KSPCG,           //Solver type
      PCGAMG,          //Preconditioner type
      1.0e-6,          //Relative residual tolerance
      1000);            //Max iterations

  // Solve
  KSPSolve(petsc_solver.ksp,b,x);

  chi::log.Log() << "Done solving";

  // Extract PETSc vector
  std::vector<double> field(num_local_dofs, 0.0);
  sdm.LocalizePETScVector(x,field,OneDofPerNode);

  // Clean up
  KSPDestroy(&petsc_solver.ksp);

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);

  chi::log.Log() << "Done cleanup";

  // Create Field Function
  auto ff = std::make_shared<chi_physics::FieldFunction>(
    "Phi",
    sdm_ptr,
    chi_math::Unknown(chi_math::UnknownType::SCALAR)
  );

  ff->UpdateFieldVector(field);

  ff->ExportToVTK("CodeTut1_FV");

  chi::Finalize();
  return 0;
}
\endcode
*/
