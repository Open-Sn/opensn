/**\page DevManMeshes The structure of meshes
 *
\section devman4_sec0 General structure of meshes

Mesh generation can be classified into 4 basic levels difficulty.
-# Orthogonal. The simplest mesh generation is undoubtedly the generation of orthogonal
meshes.
-# Reading meshes. Minimal repackaging of pre-generated meshes.
-# Semi-generated meshes. Extruded meshes and other similar operations.
-# Internally generated meshes. Pipe dream for now, but having the ability
   to generate tetrahedral meshes would open a lot of doors. From there
   I would like to have polyhedral meshes generation then hexahedral mesh
   generation.

\subsection devman_meshes_sec0_0 What is the grid?

The grid is an instance of the chi_mesh::MeshContinuum class. It has 3
members of primary interest namely:
- chi_mesh::MeshContinuum::vertices containing pointers to all the nodes of type
  chi_mesh::Node. The container is an indexed map using vertex global ids.
- chi_mesh::MeshContinuum::local_cells containing local cells of type
  chi_mesh::Cell. The container has an iterator available and is indexed by
  cell local index.
- chi_mesh::MeshContinuum::cells handles indexing by global cell ids. This
  also facilitates the addition of cells to the local_cell storage.

## _

\subsection devman_meshes_sec0_1 Easily access an existing grid

We will explain the different elements of a mesh in more detail in sections
below but as a quick primer we can say that mesh entities are loaded into
chi_mesh::MeshHandler.

You can get the current mesh-handler using:

 \code
    auto cur_handler = chi_mesh::GetCurrentHandler();
 \endcode

The computational grid is contained in a chi_mesh::MeshContinuum object which
you can obtain with:

\code
auto grid_ptr = cur_handler->GetGrid();
\endcode

If there is no existing grid then one can be created as detailed in 

## _
\section devman_meshes_sec_1 More detail on Mesh data structures

\subsection devman_meshes_sec1_0 chi_mesh::MeshHandler
Meshes and mesh operations are all handled under the umbrella of a
chi_mesh::MeshHandler. Mesh handlers are loaded onto the global variable
`chi_meshhandler_stack` and the "current" handler is tracked by
another global variable, `chi_current_mesh_handler`.
 If any mesh operation is performed it should push
newly created objects to a stack in a handler. The figure below shows
a generalized architecture.

\image html MeshOverview.png "Figure 1: Overview of the mesh hierarchy." width=700px

Even though we already showed how to obtain the handler, the "current"
 handler can always be obtained with:

 \code
    auto cur_handler = chi_mesh::GetCurrentHandler();
 \endcode

## _

\subsection devman_meshes_sec1_2 chi_mesh::SurfaceMesher and chi_mesh::VolumeMesher

Each mesh-handler is outfitted with one surface mesher and one volume
mesher, both initially undefined. The surface meshing step can be thought of
a preprocessing step. For those who are used to STAR-CCM+, this step is
analogous to the remeshing of a CAD surface into a format more conducive to
volume meshing.

There are various types of surface meshers:
 - chi_mesh::SurfaceMesherPredefined. Pass through pre-processor. Does not
   perform any meshing.
 - chi_mesh::SurfaceMesherDelaunay. Remeshes a surface mesh using Delaunay
   triangulation.
 - chi_mesh::SurfaceMesherTriangle. Remeshes a surface mesh using Triangle 1.6.

Similarly there are also various types of volume meshers:
 - chi_mesh::VolumeMesherLinemesh1D. Converts line meshes into slabs.
 - chi_mesh::VolumeMesherPredefined2D. Converts predefined surface meshes into
   2D triangles, quadrilaterals or polygons.
 - chi_mesh::VolumeMesherExtruder. Extrudes predefined surface meshes into
   3D triangular prisms, hexahedrals or polyhedrons.
 - chi_mesh::VolumeMesherPredefined3D. Converts loaded 3D meshes to
   3D tetrahedrals, hexahedrals or polyhedrons.
 - chi_mesh::VolumeMesherPredefinedUnpartitioned. Converts a ligthweight
   unpartitioned mesh to a proper full detail partitioned mesh. This mesher
   allows much flexibility for reading meshes from external sources.

Surface meshers and volume meshers are assigned to a handler as:

\code
cur_handler->surface_mesher = new chi_mesh::SurfaceMesherPredefined;
cur_handler->volume_mesher = new chi_mesh::VolumeMesherPredefined3D;
\endcode

To execute the surface meshers simply do:

\code
cur_handler->surface_mesher->Execute();
cur_handler->volume_mesher->Execute();
\endcode

## _

\subsection devman_meshes_sec1_3 chi_mesh::MeshContinuum (or if you like ... THE GRID!)

A chi_mesh::MeshContinuum object is the business-end of meshes. The execution
of a volume mesher ultimately results in the creation of a grid. To obtain
a reference to a grid simply execute:

\code
auto grid = cur_handler->GetGrid();
\endcode

This command will retrieve the latest grid from the latest region within the
current mesh handler.

## _

\subsection devman_meshes_sec1_4 chi_mesh::Cell

Cells in Chi-Tech are the basic building blocks for mesh-based scientific
computing. Some of the mesh types are shown in Figure 2 below and are defined
by the enumeration they hold as defined by chi_mesh::CellType. The cell types
supported right now are:
 - chi_mesh::CellType::SLAB
 - chi_mesh::CellType::TRIANGLE
 - chi_mesh::CellType::QUADRILATERAL
 - chi_mesh::CellType::POLYGON
 - chi_mesh::CellType::TETRAHEDRON
 - chi_mesh::CellType::HEXAHEDRON
 - chi_mesh::CellType::POLYHEDRON


\image html TypesOfCells.png "Figure 2: Types of cells." width=500px

## _

\subsection devman_meshes_sec1_5 Accessing cells (the pain of parallel programs)

Cells live in the `local_cells` member of a grid, which is of object type
chi_mesh::MeshContinuum::LocalCells, under the auspices of either
`native_cells` or `foreign_cells`. Local cells are guaranteed to be fully
defined (i.e. not ghost cells) whilst non-local cells are most likely to be
ghost-cells. The fastest way to access a cell is by means of its local id

\code
auto cell = grid->local_cells[cell_local_id];
\endcode

The `local_cells` object also facilitates an iterator.

\code
for (auto cell = grid->local_cells.begin();
     cell != grid->local-cells.end();
     ++cell)
{//do stuff}
\endcode

and also a range based iterator

\code
for (auto& cell : grid->local_cells)
{ //do stuff}
\endcode

Alternatively, a more expensive way to access a cell is by means of its
global index via the grid's `cells` member. This member is a utility object
that will search the `native_cells` and `foreign_cells` for the cell with
the associated global id

\code
auto cell = grid->cells[cell_global_id];
\endcode

Because of the search operation it is not recommended to repeatedly access
cells by their global ids.

\image html CellMapping.png "Figure 3: Cell mapping logic." width=600px
 */