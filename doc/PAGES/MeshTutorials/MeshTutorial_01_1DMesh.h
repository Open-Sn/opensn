/** \page MeshTutorial_01 Mesh Tutorial 1: Simple 1D Meshes

### It doesn't get simpler than this

Simply define the nodes of the 1D mesh in a lua table as shown below and call
the function chiMeshCreate1DSlabMesh(). As per usual, a mesh-handler needs to be
created through chiMeshHandlerCreate() and the whole process is concluded with
a call to chiVolumeMesherExecute().

\code
chiMeshHandlerCreate()
nodes={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned1DOrthoMesh(nodes)
\endcode

*/
