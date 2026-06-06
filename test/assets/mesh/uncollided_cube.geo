SetFactory("OpenCASCADE");

If (!Exists(mesh_size))
  mesh_size = 0.0064;
EndIf

Box(1) = {0.0, 0.0, 0.0, 0.032, 0.032, 0.032};
Physical Volume("domain", 1) = {1};

Mesh.MeshSizeMin = mesh_size;
Mesh.MeshSizeMax = mesh_size;
Mesh.Algorithm3D = 1;
