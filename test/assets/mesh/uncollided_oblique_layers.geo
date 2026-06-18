SetFactory("Built-in");

h = 0.06;
Point(1) = {-1.0, -1.0, 0.0, h};
Point(2) = {0.1, -1.0, 0.0, h};
Point(3) = {0.6, -1.0, 0.0, h};
Point(4) = {1.0, -1.0, 0.0, h};
Point(5) = {1.0, 1.0, 0.0, h};
Point(6) = {-0.1, 1.0, 0.0, h};
Point(7) = {-0.6, 1.0, 0.0, h};
Point(8) = {-1.0, 1.0, 0.0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 7};
Line(10) = {3, 6};

Curve Loop(1) = {1, 9, 7, 8};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 10, 6, -9};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 4, 5, -10};
Plane Surface(3) = {3};

Physical Surface("material_0", 1) = {1};
Physical Surface("material_1", 2) = {2};
Physical Surface("material_2", 3) = {3};
