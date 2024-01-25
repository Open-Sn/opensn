# Gallery

## 1. Capable of sweeps on Polyhedral meshes
Sphere embedded within a box:
- Concave cells induce cyclic dependencies between cells.
- Non-cutting Volumetric partitioning induces cyclic dependencies
  between processors.
- Data structures allows for sweeping will all forms of cycles.

![yes](images/CoolPics/SOD_threshold.png)
3D Polyhedral mesh generated with STAR-CCM+.
![yes](images/CoolPics/SOD_slice.png)
Slice of the solution.
![yes](images/CoolPics/SOD_partitioning.png)
Arbitrary, non-cutting, KBA-style partitioning used.

## 2. C5G7 Criticality Benchmark with 768 processors
Very old results (circa 2020) for the reactor benchmark C5G7.
- these results did not use acceleration techniques
- **To be updated in Spring 2024**


- 7 energy groups
- 200 directions
- 454,491 cells
- Ran on 768 processors
- Took only 18 minutes to complete (**less than 2 min today**)
- Used 584 GB of memory
- `k_eff` within 100 pcm (**much less today**)


![yes](images/CoolPics/C5G7_materials.png)
Closeup view of the mesh used. Colors represent materials.
![yes](images/CoolPics/C5G7_group0.png)
Energy group 0 solution.
![yes](images/CoolPics/C5G7_group6.png)
Energy group 6 solution.
![yes](images/CoolPics/C5G7_partition768.png)
ParMETIS partitioning of the mesh (768 processors).
![yes](images/CoolPics/C5G7_partition768b.png)
Closeup of the ParMETIS partitioning with the mesh visible.

## 2. Real world simulations
Center for Exascale Radiation Transport (CERT) simulated, and
compared to experiment, a graphite pile with a high energy neutron
source. This simulation used:
- ~172 energy groups.
- over 3000 directions.
- ~500k cells.
- Over 100k processors for some simulations.

![yes](images/CoolPics/CERTSim.png)
