/** \page MeshTutorial_03 Mesh Tutorial 3: A 3D Orthogonal Mesh

### Again ... it doesn't get simpler than this

Simply define the x,-y and z-nodes each in a lua table as shown below and call
the function chiMeshCreate3DOrthoMesh(). As per usual, a mesh-handler needs to be
created through MeshHandlerCreate() and the whole process is concluded with
a call to chiVolumeMesherExecute().

\code
MeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesz={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
MeshCreateUnpartitioned3DOrthoMesh(nodesx,nodesy,nodesz)
chiVolumeMesherExecute();
\endcode

*/
