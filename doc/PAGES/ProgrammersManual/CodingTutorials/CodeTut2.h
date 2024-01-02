/** \page CodeTut2 Coding Tutorial 2 - Using ghosted vectors

## Table of contents
- \ref CodeTut2Sec1
- \ref CodeTut2Sec2
- \ref CodeTut2Sec3
- \ref CodeTut2Sec4
- \ref CodeTut2Sec5
- \ref CodeTut2Sec6

\section CodeTut2Sec1 1 Ghost cell's in ChiTech
For the discussions that follow refer to the schematic in the figure below.
 \image html CodingTutorials/GhostCells.png width=800px

The concept of ghost cells is quite important in ChiTech, as is the case with
any scientific computing library that runs in parallel. Ghost cells arises
because it is often required to partition a mesh among different processors in
order to more effectively leverage the distributed memory models of HPCs.

<b>Definition: In ChiTech a ghost cell is a cell that shares a face or a vertex with cells
partitioned to be local to a given processor</b>. We can observe these cells in the
schematic above. This definition requires us to store the first layer of cells beyond
the locally owned cells, as well as their accompanying vertices. This allows us
to have the complete geometrical information of those ghost cells. Because some numerical
schemes often require us to map indices to ghost nodes/cells, this storage scheme,
is considered very much necessary.

In this coding tutorial we will be looking at the basics of using the spatial
discretization methods (SD methods) to gain access to ghost data, which can be
very difficult.

\section CodeTut2Sec2 2 Computing the gradient of a Finite Volume Scalar
In the finite volume SD scheme the computation of the gradient of a scalar quantity
can often be required. We will use \ref CodeTut1 as a reference for this tutorial
and assume that we are adding code to this tutorial.

Given that we want
\f$
\boldsymbol{\nabla} \phi
\f$,
how do we do this in ChiTech for the Finite Volume SD? Well the first thing to do
is to note that, from Gauss's Divergence Theorem,
\f[
\int_V \boldsymbol{\nabla} \phi dV = \int_S \mathbf{n} \phi dA
\f]
Now, if the storage of \f$ \boldsymbol{\nabla} \phi \f$ also follows the Finite
Volume SD, then \f$ \boldsymbol{\nabla} \phi \f$ will be constant on cell \f$ P \f$,
therefore
\f[
\int_S \mathbf{n} \phi dA = V_P \biggr( \boldsymbol{\nabla} \phi \biggr)_P
\f]
from which we can discretize the surface integral to get
\f[
\biggr( \boldsymbol{\nabla} \phi \biggr)_P = \frac{1}{V_P}
\sum_f A_f \mathbf{n}_f \phi_f
\f]
Now our first problem is to find \f$ \phi_f \f$, considering additionally that
we might not have an equally spaced mesh. Consider the diagram below

 \image html CodingTutorials/SqueezedHexes.png width=150px

we can now define \f$ \phi_f \f$ as a geometrically weighted average of
\f$ \phi_N \f$ and \f$ \phi_P \f$ by using the face centroid, \f$ \mathbf{x}_f \f$.
Therefore,
\f[
\phi_f = \frac{1}{|| \mathbf{x}_N - \mathbf{x}_P ||}
\biggr(
(|| \mathbf{x}_N - \mathbf{x}_f ||) \phi_P +
(|| \mathbf{x}_f - \mathbf{x}_P ||) \phi_N
\biggr)
\f]

Here we can see that we would have a problem if \f$ \phi_N \f$ is not locally
available, i.e., it is a ghost value. Therefore we need a means to obtain this
value.




\section CodeTut2Sec3 3 Making a ghosted vector and communicating ghosted data
Ghost value communication is handled via a `chi_math::VectorGhostCommunicator`
object. These objects require quite a bit of MPI communication during initializing
but thereafter can be re-used efficiently on different vectors (as long as they
share the intended structure of the ghost-communicator).

In order to gain access to `chi_math::VectorGhostCommunicator` we first include
the header
\code
#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"
\endcode

Next we add the following code:
\code
// Make ghosted vectors
std::vector<int64_t> ghost_ids = sdm.GetGhostDOFIndices(OneDofPerNode);

chi_math::VectorGhostCommunicator vgc(num_local_dofs,
                                      num_globl_dofs,
                                      ghost_ids,
                                      mpi.comm);
std::vector<double> field_wg = vgc.MakeGhostedVector(field);

vgc.CommunicateGhostEntries(field_wg);
\endcode
The first item on the agenda is "what ghost data to include in a local vector?".
To obtain the relevant ghost DOF-ids we ask the SD to provide us with the necessary
information with a call to `chi_math::spatial_discretization::GetGhostDOFIndices`.
Note this routine requires as a parameter a structure of unknowns, which is why
we supply the configuration `OneDofPerNode`.

Next we construct the `chi_math::VectorGhostCommunicator`. The constructor of
this object requires the number of local DOFs, the number of global DOFs, and
the global-ids of the ghosted data. Finally a communicator specification is
required. The construction of this object initializes data structures that will
make repeated communications very efficient.

Given that we now have the ghost-communicator we can now create a vector that
has enough space for the local data AND the ghost data. This is done with a call
to `chi_math::VectorGhostCommunicator::MakeGhostedVector` which is overloaded
to accept a vector containing only local data. This allows us, in this case,
to make a ghosted vector, `field_wg`, with the local entries the same as `field` but
with additional zeros representing the ghost values.

Finally we communicate ghost entries, of a given vector, with a call to
`chi_math::VectorGhostCommunicator::CommunicateGhostEntries`. Next we will look
at how we are going to access the ghost values.





\section CodeTut2Sec4 4 Accessing ghosted data
As a reminder, consider that a call to `chi_math::spatial_discretization::MapDOF`
provides the <b>global-id</b> of a DOF whereas
`chi_math::spatial_discretization::MapDOFLocal` provides the <b>local-id</b> of a DOF.
Ghost data is, by its nature, not normally a local quantity. Therefore, one would
assume that a call to `MapDOFLocal` would provide an invalid address for
ghost-DOFs, however, there is a catch. The id provided by `MapDOFLocal` maps to a
location aligned with the call to `GetGhostDOFIndices`. For example, consider there
are 100 local-DOFs, and 15 ghost-DOFs on a given processor. A call to
`GetGhostDOFIndices` would return 15 indices. Now lets say we want the local-id of
the 9-th ghost-DOF. In this case a call to `MapDOFLocal` will return 100+8, because
the 9-th ghost-DOF has an index of 8 in the list of ghost-ids.

We will now use this concept in the computation of the gradient of \f$ \phi \f$.

First we create a vector that can store the x, y and z components of the gradient.
The size of this vector, as well as mappings into it, will be provided by the SD
and therefore we need to make a `chi_math::UnknownManager` to dictate the structure
of the unknowns. Therefore we do the following:
\code
chi_math::UnknownManager grad_uk_man(
    {chi_math::Unknown(chi_math::UnknownType::VECTOR_3)});
\endcode

Next we get the size of the vector and instantiate it
\code
const size_t num_grad_dofs = sdm.GetNumLocalDOFs(grad_uk_man);

std::vector<double> grad_phi(num_grad_dofs, 0.0);
\endcode

Since we already have ghosted data we are now ready to compute the gradient on
a cell-by-cell basis. We do this with the following code:
\code
for (const auto& cell : grid.local_cells)
{
  const auto& cell_mapping = sdm.GetCellMapping(cell);

  const int64_t imap = sdm.MapDOFLocal(cell, 0);
  const double  phi_P = field_wg[imap];

  const auto& xp = cell.centroid;

  auto grad_phi_P = chi_mesh::Vector3(0,0,0);

  size_t f=0;
  for (const auto& face : cell.faces)
  {
    const auto& xf = face.centroid;
    const auto  Af = cell_mapping.FaceArea(f) * face.normal;

    double phi_N = 0.0;
    auto xn = xp + 2*(xf-xp);

    if (face.has_neighbor)
    {
      const auto& adj_cell = grid.cells[face.neighbor_id];
      const int64_t nmap = sdm.MapDOFLocal(adj_cell, 0);
      phi_N = field_wg[nmap];

      xn = adj_cell.centroid;
    }

    grad_phi_P += Af * ((xn-xf).Norm()*phi_P + (xf-xp).Norm()*phi_N)/
                  (xn-xp).Norm();
    ++f;
  }//for face
  grad_phi_P /= cell_mapping.CellVolume();

  const int64_t xmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 0);
  const int64_t ymap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 1);
  const int64_t zmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 2);

  grad_phi[xmap] = grad_phi_P.x;
  grad_phi[ymap] = grad_phi_P.y;
  grad_phi[zmap] = grad_phi_P.z;
}//for cell
\endcode
As per usual we loop over the local cells. Once we have a cell reference we
obtain a cell-mapping. We then immediately grab the present cell's \f$ \phi_P \f$
value by first mapping its address (reminder, the second paramater of `MapDOFLocal`
is the node-number, which is zero for the Finite Volume SD), then accessing the
ghosted field vector, `field_wg`.

We then initialize an, initially zero, three-component vector representing the
gradient for cell \f$ P \f$. Next we proceed to loop over the faces of the cell
where we obtain the neighboring cell's \f$ \phi_N \f$ value (provided that
`face.has_neighbor`, otherwise defaulting to the boundary condition of zero)
with the code
\code
const auto& adj_cell = grid.cells[face.neighbor_id];
const int64_t nmap = sdm.MapDOFLocal(adj_cell, 0);
phi_N = field_wg[nmap];
\endcode
Notice here the use of `MapDOFLocal`. This only works because `field_wg` will have
a value for the case when `MapDOFLocal` is asked to map a ghosted DOF.

Finally we compute the gradient from the formula
\f[
\phi_f = \frac{1}{|| \mathbf{x}_N - \mathbf{x}_P ||}
\biggr(
(|| \mathbf{x}_N - \mathbf{x}_f ||) \phi_P +
(|| \mathbf{x}_f - \mathbf{x}_P ||) \phi_N
\biggr)
\f]
with the code
\code
grad_phi_P += Af * ((xn-xf).Norm()*phi_P + (xf-xp).Norm()*phi_N)/
                    (xn-xp).Norm();
\endcode
which is divided by the volume at the termination of the faces-loop.

The final step for a cell is to appropriately store this gradient, therefore
we map each component of the gradient and appropriately set the `grad_phi`
vector,
\code
const int64_t xmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 0);
const int64_t ymap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 1);
const int64_t zmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 2);

grad_phi[xmap] = grad_phi_P.x;
grad_phi[ymap] = grad_phi_P.y;
grad_phi[zmap] = grad_phi_P.z;
\endcode





\section CodeTut2Sec5 5 Visualizing the gradient
Visualizing the gradient is conceptually simple. We simply create a
field-function that will export all three the components of the gradient together.
Then we use paraview to make arrow glyphs for visualization.
\code
// Create Field Function
auto ff_grad = std::make_shared<chi_physics::FieldFunction>(
  "GradPhi",
  sdm_ptr,
  chi_math::Unknown(chi_math::UnknownType::VECTOR_3)
);

ff_grad->UpdateFieldVector(grad_phi);

ff_grad->ExportToVTK("CodeTut2_FV_grad");
\endcode

Within paraview we can use the `Calculator` filter to construct vector quantities
from the 3 components after which we can apply a `Glyph` filter to obtain a
visualization as shown below

\image html CodingTutorials/Tut2Gradient.png width=700px






\section CodeTut2Sec6 The complete program

\code
#include "chi_runtime.h"
#include "chi_log.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/FiniteVolume/fv.h"
#include "math/PETScUtils/petsc_utils.h"

#include "physics/FieldFunction/fieldfunction2.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"

int main(int argc, char* argv[])
{
  chi::Initialize(argc,argv);
  chi::RunBatch(argc, argv);

  chi::log.Log() << "Coding Tutorial 2";

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

  ff->ExportToVTK("CodeTut2_FV");

  // Make ghosted vectors
  std::vector<int64_t> ghost_ids = sdm.GetGhostDOFIndices(OneDofPerNode);

  chi_math::VectorGhostCommunicator vgc(num_local_dofs,
                                        num_globl_dofs,
                                        ghost_ids,
                                        mpi.comm);
  std::vector<double> field_wg = vgc.MakeGhostedVector(field);

  vgc.CommunicateGhostEntries(field_wg);

  // Setup gradient unknown structure
  chi_math::UnknownManager grad_uk_man(
    {chi_math::Unknown(chi_math::UnknownType::VECTOR_3)});

  const size_t num_grad_dofs = sdm.GetNumLocalDOFs(grad_uk_man);

  std::vector<double> grad_phi(num_grad_dofs, 0.0);

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);

    const int64_t imap = sdm.MapDOFLocal(cell, 0);
    const double  phi_P = field_wg[imap];

    const auto& xp = cell.centroid;

    auto grad_phi_P = chi_mesh::Vector3(0,0,0);

    size_t f=0;
    for (const auto& face : cell.faces)
    {
      const auto& xf = face.centroid;
      const auto  Af = cell_mapping.FaceArea(f) * face.normal;

      double phi_N = 0.0;
      auto xn = xp + 2*(xf-xp);

      if (face.has_neighbor)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id];
        const int64_t nmap = sdm.MapDOFLocal(adj_cell, 0);
        phi_N = field_wg[nmap];

        xn = adj_cell.centroid;
      }

      grad_phi_P += Af * ((xn-xf).Norm()*phi_P + (xf-xp).Norm()*phi_N)/
                    (xn-xp).Norm();
      ++f;
    }//for face
    grad_phi_P /= cell_mapping.CellVolume();

    const int64_t xmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 0);
    const int64_t ymap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 1);
    const int64_t zmap = sdm.MapDOFLocal(cell, 0, grad_uk_man, 0, 2);

    grad_phi[xmap] = grad_phi_P.x;
    grad_phi[ymap] = grad_phi_P.y;
    grad_phi[zmap] = grad_phi_P.z;
  }//for cell

  // Create Field Function
  auto ff_grad = std::make_shared<chi_physics::FieldFunction>(
    "GradPhi",
    sdm_ptr,
    chi_math::Unknown(chi_math::UnknownType::VECTOR_3)
  );

  ff_grad->UpdateFieldVector(grad_phi);

  ff_grad->ExportToVTK("CodeTut2_FV_grad");

  chi::Finalize();
  return 0;
}
\endcode

*/
