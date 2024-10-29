## OpenSn Unstructured Mesh Node-to-Node Weak Scaling Study

This is a node-to-node weak scaling study for OpenSn with the following
characteristics:

- Unstructured mesh
- 64 ranks per node
- 256 cells per core
- 448 angles
- 64 energy groups,
- PWLD spatial discretization

The `.geo` file needs to be processed by `Gmsh` to create the necessary `.msh`
file. 

`slurm.sh` is provided as an example job script. For each job in the study, 
both the number of nodes and the name of the input file need to be modifed
in the job script.
