// Parameters
cube_size = 10.0;          // Side length of the cube

// Calculate characteristic length based on the desired number of tetrahedra
approx_element_size = cube_size / 15;

// Define cube corner points
Point(1) = {0, 0, 0, approx_element_size};
Point(2) = {cube_size, 0, 0, approx_element_size};
Point(3) = {cube_size, cube_size, 0, approx_element_size};
Point(4) = {0, cube_size, 0, approx_element_size};
Point(5) = {0, 0, cube_size, approx_element_size};
Point(6) = {cube_size, 0, cube_size, approx_element_size};
Point(7) = {cube_size, cube_size, cube_size, approx_element_size};
Point(8) = {0, cube_size, cube_size, approx_element_size};

// Define the edges of the cube (lines)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

// Define surface loops for each face of the cube
Line Loop(1) = {1, 2, 3, 4};           // Bottom face
Line Loop(2) = {1, 6, -9, -5};         // Side face between points 1, 2, 6, 5
Line Loop(3) = {2, 7, -10, -6};        // Side face between points 2, 3, 7, 6
Line Loop(4) = {3, 8, -11, -7};        // Side face between points 3, 4, 8, 7
Line Loop(5) = {4, 5, -12, -8};        // Side face between points 4, 1, 5, 8
Line Loop(6) = {9, 10, 11, 12};        // Top face

// Create plane surfaces from the surface loops
Plane Surface(1) = {1};  // Bottom surface
Plane Surface(2) = {2};  // Side surface 1
Plane Surface(3) = {3};  // Side surface 2
Plane Surface(4) = {4};  // Side surface 3
Plane Surface(5) = {5};  // Side surface 4
Plane Surface(6) = {6};  // Top surface

// Create a surface loop to define the boundary of the volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};

// Create the volume using the surface loop
Volume(1) = {1};

// Mesh settings
Mesh.CharacteristicLengthMin = approx_element_size;
Mesh.CharacteristicLengthMax = approx_element_size;
Mesh.ElementOrder = 1;    // First-order tetrahedra
Mesh 3;                   // Generate the 3D tetrahedral mesh

// Save and partition the mesh
Physical Volume("CubeVolume") = {1};
Physical Surface("BottomFace") = {1};
Physical Surface("TopFace") = {6};

