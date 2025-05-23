## OpenSn Unstructured Mesh Node-to-Node Strong Scaling Study

This is a node-to-node strong scaling study for OpenSn with the following
characteristics:

- Unstructured mesh
- 268,938 cells (~4,200 cells/rank at 1 node to ~8 cells/rank at 512 nodes)
- 64 ranks per node
- 448 angles
- 64 energy groups,
- PWLD spatial discretization

The Python file `generate_scaling_study.py` can be used to generate the mesh
files, inputs, and SLURM job scripts for the scaling study. Gmsh must be 
installed to generate the mesh files. Various configuration options including
the location of the Gmsh binary, the location of the OpenSn binary, tasks per
node, number of nodes, etc. are located at the bottom of the 
`generate_scaling_study.py` script.
