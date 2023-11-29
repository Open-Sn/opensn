/** \page MeshTutorial_02 Mesh Tutorial 2: A 2D Orthogonal Mesh

### Again ... it doesn't get simpler than this

Simply define the x- and y-nodes each in a lua table as shown below and call
the function chiMeshCreate2DOrthoMesh(). As per usual, a mesh-handler needs to be
created through MeshHandlerCreate() and the whole process is concluded with
a call to chiVolumeMesherExecute().

\code
MeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
MeshCreateUnpartitioned2DOrthoMesh(nodesx,nodesy)
chiVolumeMesherExecute();
\endcode

*/
