/**\page DevManUKMan The chi_math::UnknownManager
## Table of contents
- \ref devmanA3_sec0
- \ref devmanA3_sec1
- \ref devmanA3_sec2
- \ref devmanA3_sec3
- \ref devmanA3_sec4


\section devmanA3_sec0 Introduction
As a developer you might come across code such as
\code
imap[i] = sdm.MapDOF(cell, i);
\endcode
or its more elaborate cousin
\code
const int64_t imap = m_sdm.MapDOF(cell,i,m_uk_man,0,g);
\endcode

The two function calls, `MapDOF` with 2 arguments and `MapDOF` with 5 arguments,
are both methods of the base-class `chi_math::SpatialDiscretization`.

The two overloaded methods, `chi_math::spatial_discretization::MapDOF` return the
global index of a degree-of-freedom (DOF). The counterparts to `MapDOF` are the
`MapDOFLocal` methods that provide the local index of a DOF.

These functions are called with either 2 arguments
- A cell reference
- A cell-node local-index

or with 5 arguments
- A cell reference
- A cell-node local-index
- An unknown-manager (`chi_math::UnknownManager`)
- An unknown id
- A component id

\section devmanA3_sec1 1 Cells vs Vertices vs Nodes vs Degree-of-Freedom (DOF)
These comparisons are conceptually simple. The geometry of a cell is described
by its vertices (see Figure 1 below). The mesh cannot be completely defined
without the vertices, therefore, cells always go hand-in-hand with vertices.

A spatial discretization (e.g., Finite Volume, PWL, Lagrange Q1, Lagrange Q2)
places nodes on the mesh that may or may not coincide with some of the vertices.
For example, the Finite Volume spatial discretization places a single node for
each cell, at the cell's centroid (hence no coincidence with any vertex), whilst
a Lagrange Q2 spatial discretization has some of its nodes on the vertices while
other are not (see Figure 1 below).

\image html UkManVertsNodes.png "Figure 1: Cells, vertices and nodes." width=600px

Where degrees-of-freedom differ from nodes is that each node can be stacked with
a number of unknowns, e.g., one might have to solve a number of physical
variables in a multiphysics simulation involving 2D velocity, \f$ u_x \f$ and
\f$ u_y \f$, a three energy group flux, \f$ \phi_0, \ \phi_1, \ \phi_2 \f$, and
pressure, \f$ p \f$. One can stack these variables in any form, for example:
- Unknown 0, a 6 component vector containing
\f$ u_x, \ u_y, \ \phi_0, \ \phi_1, \ \phi_2 \f$, and \f$ p \f$

or more elegantly
- Unknown 0, a 2 component vector containing \f$ u_x \f$ and \f$ u_y \f$
- Unknown 1, a 3 component vector containing \f$ \phi_0, \ \phi_1, \ \phi_2 \f$
- Unknown 2, a single scalar (1 component) containing just \f$ p \f$

The latter essentially means that each node can have a number of unknowns stacked
onto it as shown in Figure 2 below.

\image html UkManNodeDOFs.png "Figure 2: Unknowns/DOFs stacked onto a node" width=300px

The individual components of the unknowns are unique for every node and are
called the degrees-of-freedom (DOFs).

\section devmanA3_sec2 2 One DOF per node
All the spatial discretizations have a default unknown-manager, defined during
initialization, called `UNITARY_UNKNOWN_MANAGER`. This allows the spatial
discretization to map indices as shown in Figure 3 below

\image html UkMan-OneDOFPN.png "Figure 3: Mapping One DOF per Node" width=800px

\section devmanA3_sec3 3 Multiple DOFs Nodal ordering
With multiple unknowns like the 2D velocity, multigroup flux, and pressure
example above we can create an unknown manager with this structure as follows
\code
chi_math::UnknownManager uk_man;
uk_man.AddUnknown(chi_math::UnknownType::VECTOR_2);    //u_x, u_y
uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, 3); //phi_0 to phi_2
uk_man.AddUnknown(chi_math::UnknownType::SCALAR);      //pressure
\endcode

Using this unknown manager to map DOFs is then done as shown below in Figure 4.

\image html UkMan-NODAL.png "Figure 4: Mapping with a nodal ordering" width=800px

Notice here that the ordering is per-node, i.e., for each node the DOFs all the
DOFs are stacked into a vector before moving to a new node. This is the default
ordering in ChiTech.

\section devmanA3_sec4 4 Multiple DOFs Block ordering
To create a Block-ordered vector mapping the only modification needed is in the
constructor of the unknown manager:
\code
chi_math::UnknownManager uk_man(chi_math::UnknownStorageType::BLOCK);
uk_man.AddUnknown(chi_math::UnknownType::VECTOR_2);    //u_x, u_y
uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, 3); //phi_0 to phi_2
uk_man.AddUnknown(chi_math::UnknownType::SCALAR);      //pressure
\endcode

Now that we specified that we should use `chi_math::UnknownStorageType::BLOCK`
the unknown manager maps as shown in Figure 5 below

\image html UkMan-Block.png "Figure 5: Mapping with a block ordering" width=800px

Block ordering is likely desireable if one wishes to apply block-solves.
 */
