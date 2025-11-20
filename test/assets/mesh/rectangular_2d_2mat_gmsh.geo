// Define geometry parameters
Lx = 10;  // Length of the plane along x-axis
Ly = 5;   // Length of the plane along y-axis

// Define points
Point(1) = {0, 0, 0};            // Bottom-left corner
Point(2) = {Lx, 0, 0};           // Bottom-right corner
Point(3) = {Lx, Ly/2, 0};        // Right mid-point (interface between two materials)
Point(4) = {0, Ly/2, 0};         // Left mid-point (interface between two materials)
Point(5) = {Lx, Ly, 0};          // Top-right corner
Point(6) = {0, Ly, 0};           // Top-left corner

// Define lines
Line(1) = {1, 2};                // Bottom edge
Line(2) = {2, 3};                // Right edge (bottom material)
Line(3) = {3, 4};                // Horizontal interface between the two materials
Line(4) = {4, 1};                // Left edge (bottom material)
Line(5) = {3, 5};                // Right edge (top material)
Line(6) = {5, 6};                // Top edge
Line(7) = {6, 4};                // Left edge (top material)

// Define two line loops, one for each material
Line Loop(1) = {1, 2, 3, 4};    // Bottom material (closed loop)
Line Loop(2) = {5, 6, 7, -3};   // Top material (closed loop, reuses interface)

// Define two surfaces, one for each material
Plane Surface(1) = {1};          // Surface for bottom material
Plane Surface(2) = {2};          // Surface for top material

// Define physical groups to represent two materials
Physical Surface("Material 1") = {1};  // Bottom material
Physical Surface("Material 2") = {2};  // Top material

// Mesh settings (optional)
Mesh.CharacteristicLengthMin = 0.2;    // Minimum element size
Mesh.CharacteristicLengthMax = 0.5;    // Maximum element size

// Generate mesh
Mesh 2;

