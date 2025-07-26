Gallery
=======

1. Capable of sweeps on Polyhedral meshes
-----------------------------------------

Sphere embedded within a box:

- Concave cells induce cyclic dependencies between cells.
- Non-cutting volumetric partitioning induces cyclic dependencies
  between processors.
- Data structures allow sweeping with all forms of cycles.

.. figure:: images/CoolPics/SOD_threshold.png
   :alt: Polyhedral mesh threshold

   3-D polyhedral mesh generated with **STAR-CCM+**.

.. figure:: images/CoolPics/SOD_slice.png
   :alt: Solution slice

   Slice of the solution.

.. figure:: images/CoolPics/SOD_partitioning.png
   :alt: KBA-style partitioning

   Arbitrary, non-cutting, KBA-style partitioning used.


2. C5G7 Criticality Benchmark with 768 processors
-------------------------------------------------

Very old results (circa 2020) for the reactor benchmark **C5G7**.

- These results did **not** use acceleration techniques  
  (**to be updated in Spring 2024**)
- 7 energy groups
- 200 directions
- 454 491 cells
- Ran on 768 processors
- Took only 18 minutes to complete (**less than 2 min today**)
- Used 584 GB of memory
- ``k_eff`` within 100 pcm (**much less today**)

.. figure:: images/CoolPics/C5G7_materials.png
   :alt: C5G7 mesh materials

   Close-up view of the mesh used (colors represent materials).

.. figure:: images/CoolPics/C5G7_group0.png
   :alt: Energy group 0 solution

   Energy group 0 solution.

.. figure:: images/CoolPics/C5G7_group6.png
   :alt: Energy group 6 solution

   Energy group 6 solution.

.. figure:: images/CoolPics/C5G7_partition768.png
   :alt: 768-way ParMETIS partition

   ParMETIS partitioning of the mesh (768 processors).

.. figure:: images/CoolPics/C5G7_partition768b.png
   :alt: ParMETIS partition close-up

   Close-up of the ParMETIS partitioning with the mesh visible.


3. Real world simulations
-------------------------

The **Center for Exascale Radiation Transport (CERT)** simulated—and
compared to experiment—a graphite pile with a high-energy neutron source.
This simulation used:

- ~172 energy groups
- Over 3000 directions
- ~500 k cells
- Over 100 k processors for some simulations

.. figure:: images/CoolPics/CERTSim.png
   :alt: CERT simulation

   CERT simulation versus experiment.
