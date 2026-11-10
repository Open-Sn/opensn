# Kobayashi Problem 3 pMesh mesh

This OpenFOAM `polyMesh` represents the same OECD/NEA Kobayashi Problem 3
geometry and material numbering as the tetrahedral meshes in the parent
directory:

- block 1, `source`: 0--10 cm in x, y, and z;
- block 2, `duct`: the four connected 10 cm-wide dog-leg segments; and
- block 3, `shield`: the remainder of the 60 x 100 x 60 cm box.

The mesh contains 9,024 arbitrary polyhedral cells and 38,281 vertices. The
source, duct, and shield contain 35, 453, and 8,536 cells, respectively. Their
material volumes are 1,000, 15,000, and 344,000 cubic centimetres. Polygonal
faces are retained rather than subdivided into triangles. The six boundary
patches are named `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, and `zmax`.

The `unused` cell zone occupies OpenFOAM zone index zero so that the physical
zones retain the benchmark block IDs 1, 2, and 3 when OpenSn reads them. Its
single cell is also present in its physical zone, which is written later and
therefore supplies the cell's effective block ID.

## Provenance and regeneration

The base cells were generated with the GPLv3-or-later open-source
[cfMesh 1.2.0](https://sourceforge.net/projects/cfmesh/) `pMesh` executable.
The source archive used had SHA-256
`429f3fb2d73c2a9831f49a5bde1a26510847a80973f86a3b3fdac3e4101e9577`.
Gmsh 4.13.1 generated the input STL. The input surface is deliberately padded
beyond the benchmark domain, and `pMesh` stops after `templateGeneration`.
This retains its convex Voronoi cells and avoids the non-star-shaped cells that
its later smoothing stages can introduce.

From an empty working directory containing `system/meshDict` and
`system/controlDict`, generate the base mesh with:

```sh
gmsh -2 kobayashi_box.geo -format stl -o kobayashi_box.stl
pMesh
python postprocess.py constant/polyMesh final/constant/polyMesh
```

The files needed for those commands are in `generation/`. The postprocessor:

1. normalizes cfMesh's padded `-1` neighbour entries;
2. clips the Voronoi cells to the exact outer and material-interface planes;
3. resolves one-to-many faces created by clipping;
4. retains consistently oriented polygonal faces; and
5. writes OpenFOAM boundary patches and cell zones in the form read natively
   by OpenSn.

The checked-in mesh is covered by the
`uncollided_3d_polyhedral_reflecting.py` regression. OpenSn reads it natively,
checks the three block volumes against the exact benchmark volumes, and
compares a reflected point-source calculation with the analytic image-source
solution. An eight-thread uncollided Kobayashi solve, including reflected image
sources on all three symmetry planes, also completed successfully on this mesh.
