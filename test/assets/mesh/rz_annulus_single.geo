// RZ case 04: annulus-style rectangle r0<=r<=r1, 0<=z<=2
r0 = 0.25;
r1 = 1.0;
z0 = 0.0;
z1 = 2.0;
lc = 0.10;

Point(1) = {r0, z0, 0, lc};
Point(2) = {r1, z0, 0, lc};
Point(3) = {r1, z1, 0, lc};
Point(4) = {r0, z1, 0, lc};

Line(1) = {1, 2}; // zmin
Line(2) = {2, 3}; // rmax
Line(3) = {3, 4}; // zmax
Line(4) = {4, 1}; // rmin (non-axis)

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface(1) = {1};
Physical Curve("zmin", 101) = {1};
Physical Curve("rmax", 102) = {2};
Physical Curve("zmax", 103) = {3};
Physical Curve("rmin", 104) = {4};

Mesh 2;
