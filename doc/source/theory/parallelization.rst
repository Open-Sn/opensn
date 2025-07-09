Parallel Sweeps
===============

Introduction
------------

A parallel sweep algorithm is comprised of the following three
components:

-  **Aggregation** - The grouping of cells, angles, and groups into
   tasks. Note that in OpenSn, due to the arbitrary nature of the grid,
   we don’t aggregate cells.

-  **Partitioning** - The division of the mesh among processors.

-  **Scheduling** - The execution of available tasks in a specified
   order by each processor.

Sweep performance is profoundly affected by the choices made for
aggregation, partitioning, and scheduling. The interplay of the three
components is complex and nuanced, and the reader is referred to the
recommended references for more detailed discussions,
:cite:t:`mathis2000general`, :cite:t:`KBA`, :cite:t:`baker1998s`, :cite:t:`bailey2009analysis`, :cite:t:`pautz2002algorithm`, :cite:t:`zeyao2004parallel`, :cite:t:`hanebutte1992massively`, :cite:t:`brown1999performing`, :cite:t:`brown2000performing`, :cite:t:`adams2013provably`, :cite:t:`eff_sweeps`, :cite:t:`vermaak2021massively`.
The content below is only intended to be a very brief introduction to
the basic concepts.

Aggregation
-----------

In the previous definitions of the multigroup iterative algorithms, the
description of the solution techniques were given on a group-by-group
basis. Likewise, the solution of a single-group problem problem was
given direction-by-direction. At the level of a given group :math:`g`,
for a given direction :math:`d`, a full-domain transport sweep is
performed where the mesh is swept cell-by-cell, solving for the angular
flux :math:`\psi^g_{d,K}` in each cell :math:`K`. A prototypical
algorithm following those descriptions is given in Algorithm 1 and shows
one potential ordering of the three phase-space loops. Each loop
ordering affects code performance in different ways. For example, having
the group loop as the outermost loop has the advantage that only one
pass through the loop is needed for problems with downscattering only.
However, it also leads to a significant number of repeated geometry
calculations when sweeping an arbitrary grid. If we swap the group and
cell loops so that the group loop is now the innermost loop, we no
longer get immediate convergence in problems with downscattering only,
but we can re-use geometry information for each group and cell solve.
The optimal ordering of loops is problem-dependent, and, while we could
introduce separate compute kernels for each ordering, a better solution
is to extend our sweep algorithm to support sets of groups and angles.

.. figure:: images/algo1.png
   :scale: 50%
   :align: center


Group-sets
~~~~~~~~~~

A group-set is a collection of one or more groups and gives us the
flexibility to support different loop orderings without requiring
multiple compute kernels. Algorithm 2 extends Algorithm 1 to support
group-sets.

.. figure:: images/algo2.png
   :scale: 50%
   :align: center


If each group is its own group-set, our solution algorithm is similar to
Algorithm 1 and is optimal for problems with downscattering only. If all
of the groups are contained in a single group-set, the group loop is now
the innermost loop, and we can solve all of the groups for each cell
concurrently. This greatly increases our flop count and amortizes the
cost of querying upwind flux values and other cell-specific quantities.

Angle-sets
~~~~~~~~~~

Similar to group-sets, angle-sets allow us to group angles together that
share similar sweep orderings (they have similar task-directed graphs).
With angle-sets, it is possible to solve for multiple angular-flux
quantities concurrently on each spatial cell. Angle-sets also allow us
to cache common geometric quantities for cell solves and to optimize our
sweep-plane data structures for cache efficiency. Algorithm 3 extends
our group-set algorithm to include angle-sets.

.. figure:: images/algo3.png
   :scale: 50%
   :align: center


Partitioning and Scheduling
---------------------------

Partitioning is the process of dividing the spatial domain among
processes. For example, consider the simple 2D mesh shown in Figure 1.
This 4x4 mesh has been decomposed among four processors, with each
processor being assigned a column of four cells. This columnar
decomposition is the type of partitioning produced by the KBA algorithm.

.. figure:: images/ProcMesh.png
   :scale: 50%
   :align: center

   Partitioned, 4x4, XY, Rectangular Mesh

Now, consider the direction :math:`\vec{\Omega}` incident on the
bottom-left corner of the mesh. The sweep ordering for this quadrature
direction is represented by the task-dependency graph shown in Figure 2.
As the graph illustrates, a cell on a given level of the graph cannot be
solved until its upstream neighbors on the previous level of the graph
have been solved. The graph is processed stage by stage, with each
processor solving its “ready” cells at each stage and communicating
outgoing angular fluxes to its downstream neighbors. As a processor
finishes its tasks for one direction in an angle-set, it begins
executing its tasks for the next direction until all angles in the
angle-set have been processed. Once a processor has finished with one
angle-set/group-set combination, it can start solving the next available
angle-set/group-set. This pipelining of tasks is important to keep each
processor working for as long as possible and to maintain parallel
efficiency. In the case where multiple angle-sets can be solved
simultaneously, a scheduling algorithm is required to tell the processor
the order in which to execute the available tasks.

.. figure:: images/ProcGraph.png
   :scale: 50%
   :align: center

   Task Dependency Graph for Direction :math:`\vec{\Omega}`.


References
----------
    
.. bibliography::
   :style: unsrtalpha
   :filter: False

   KBA
   adams2013provably
   bailey2009analysis
   baker1998s
   brown1999performing
   brown2000performing
   eff_sweeps
   hanebutte1992massively
   mathis2000general
   pautz2002algorithm
   vermaak2021massively
   zeyao2004parallel
