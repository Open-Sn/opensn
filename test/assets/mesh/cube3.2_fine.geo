// Refined companion mesh to cube3.2.msh for the uncollided-flux tutorial
// (doc/source/tutorials/workflows/data_reuse/uncollided).
//
// Same 0.032 m cube domain as cube3.2.msh, graded around the tutorial's
// point-source location: 0.6 mm elements within 4 mm of the source, growing
// to 1.8 mm beyond 14 mm. The graded (rather than uniform) refinement is
// deliberate: a uniform fine mesh does not control how close the point
// source lands to a cell face, and OpenSn's PWLD point-source representation
// is sensitive to that -- an off-center or near-face placement measurably
// biases the near-source ray-traced flux (see the "point source ... lies
// only <dist> from a face" warning in uncollided_problem.cc). Grading the
// mesh around the source keeps the source-containing cell small and
// well-shaped regardless of exactly where Delaunay meshing places it.
//
// Regenerate with (gmsh 4.15.2 used originally):
//   gmsh -3 cube3.2_fine.geo -o cube3.2_fine.msh -format msh41

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 0.032, 0.032, 0.032};

source_x = 0.0102586;
source_y = 0.0114131;
source_z = 0.0146416;

p = newp;
Point(p) = {source_x, source_y, source_z, 1.0};

Field[1] = Distance;
Field[1].PointsList = {p};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.0006;
Field[2].SizeMax = 0.0018;
Field[2].DistMin = 0.004;
Field[2].DistMax = 0.014;

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MshFileVersion = 4.1;
Mesh.SaveAll = 1;
Mesh.OptimizeNetgen = 1;
