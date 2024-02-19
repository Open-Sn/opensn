# Introduction to meshing and partitioning

## Meshing

For meshing, we recommend users to utilize mesh generators that can generate files in the following formats:
- gmesh 
- vtu 
- Exodus.II
- Ensight gold
- Wavefront (obj) [only in 2D]

Meshes can be supplied with `material IDs` and `boundary IDs` to define heterogeneous
properties in the domain (i.e., materials and volumetric sources) and boundary sources. 

Support for additional mesh formats may be added in the future.

We also provide our meshing capability for orthogonal grids (in 1D, 2D, and 3D).

## Partitioning

For parallel simulations, the domain is split (partitioned) either using 
- the `Parmetis partitioner` (through PETSc) 
or 
- the `KBA partitioner` (mostly useful for orthogonal grids).
