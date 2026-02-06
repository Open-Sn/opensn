// RZ case 08: rectangular quad mesh, 0<=r<=1, 0<=z<=2
r0 = 0.0;
r1 = 1.0;
z0 = 0.0;
z1 = 2.0;

nr = 50;
nz = 100;

Point(1) = {r0, z0, 0, 1.0};
Point(2) = {r1, z0, 0, 1.0};
Point(3) = {r1, z1, 0, 1.0};
Point(4) = {r0, z1, 0, 1.0};

Line(1) = {1, 2}; // zmin
Line(2) = {2, 3}; // rmax
Line(3) = {3, 4}; // zmax
Line(4) = {4, 1}; // rmin

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = nr + 1;
Transfinite Line {2, 4} = nz + 1;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface(1) = {1};
Physical Curve("zmin", 101) = {1};
Physical Curve("rmax", 102) = {2};
Physical Curve("zmax", 103) = {3};
Physical Curve("rmin", 104) = {4};

Mesh 2;
