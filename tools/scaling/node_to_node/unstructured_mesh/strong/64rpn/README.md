## OpenSn Unstructured Mesh Node-to-Node Strong Scaling Study

This is a node-to-node strong scaling study for OpenSn with the following
characteristics:

- Unstructured mesh
- 268,938 cells (~4,200 cells/rank at 1 node to ~8 cell/rank at 512 nodes)
- 64 ranks per node
- 448 angles
- 64 energy groups,
- PWLD spatial discretization

The `.geo` file needs to be processed by `Gmsh` to create the necessary `.msh`
file. 

`slurm.sh` is provided as an example job script. For each job in
the study, only the number of nodes in the job script needs to be modifed.
