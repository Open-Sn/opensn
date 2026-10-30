================
Running Problems
================

This section explains how to run OpenSn problems once an input is written.
The examples here focus on the OpenSn console application.

Basic command-line usage
========================

An OpenSn Python input is typically run with the ``opensn`` executable and the
``-i`` option:

.. code-block:: bash

   build/python/opensn -i my_input.py

Most users run inputs this way during development and debugging. The input
script creates meshes, materials, problems, and solvers, then calls the solver
methods needed for the calculation.

.. note::

   OpenSn inputs are Python scripts. This means a problem input can contain
   loops, helper functions, conditionals, and post-processing logic in addition
   to the solver setup itself.

Running under MPI
=================

Parallel runs are launched with ``mpirun`` or the site-specific MPI launcher:

.. code-block:: bash

   mpirun -np 4 build/python/opensn -i my_input.py

The number of MPI ranks should be chosen with the mesh partitioning and problem
size in mind. Very small problems generally do not benefit from large MPI
counts, while large three-dimensional transport problems are usually intended
to run in parallel.

.. note::

   If a mesh generator or imported mesh already contains partitioning
   information, make sure that the run configuration is compatible with the
   intended decomposition.

The ``max_mpi_message_size`` problem option sets the largest MPI message, in
bytes, that an AAH sweep sends, on the CPU and on the GPU. Larger angular-flux
transfers are split into messages of at most this size. The default, 4194304
bytes (4 MiB), sends most transfers whole. Small limits split large transfers
into many messages, and each one adds posting, matching, and polling work. With
a 32 KiB limit, GPU sweeps on a few nodes ran up to 2.8 times slower. Messages
below the limit are not affected, so the default rarely needs changing.

Running with GPUs
=================

GPU acceleration is requested in the problem definition with ``use_gpus=True``.
For example:

.. code-block:: python

   phys = opensn.DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=20,
       groupsets=groupsets,
       xs_map=xs_map,
       use_gpus=True,
       sweep_type="AAH",
   )

This does not change the command-line launch pattern. A GPU-enabled input is
still run with the normal OpenSn executable, either in serial or under MPI:

.. code-block:: bash

   mpirun -np 4 build/python/opensn -i my_input.py

Important restrictions:

* OpenSn must be built with GPU support.
* Curvilinear discrete ordinates problems do not support GPU acceleration.
* Time-dependent problems do not support GPU acceleration.
* OpenSn uses one GPU per MPI rank.

Practical guidance:

* First make the problem work on the CPU path.
* Then enable ``use_gpus=True`` and keep the rest of the input as unchanged as
  possible.
* If a GPU run fails immediately, check the build configuration, and whether the
  problem is steady-state or transient.

Each MPI rank selects its GPU at startup, and PETSc uses the same GPU:

* Ranks on a node that see the same GPUs share them in turn by node-local rank.
  With four GPUs, node-local ranks 0-3 use GPUs 0-3, rank 4 uses GPU 0, and so
  on.
* A batch system or launcher that assigns GPUs to ranks does so by restricting
  the GPUs each rank sees, for example with ``CUDA_VISIBLE_DEVICES``,
  ``ROCR_VISIBLE_DEVICES``, or ``HIP_VISIBLE_DEVICES`` (Slurm's
  ``--gpu-bind``, for instance). A rank given one GPU uses it; ranks given the
  same GPUs share them in turn.
* OpenSn passes the selected GPU to PETSc as ``-device_select_cuda`` or
  ``-device_select_hip``. Do not select a different PETSc device, for example
  through ``PETSC_OPTIONS``.
* With SYCL, the SYCL runtime selects the GPU. Restrict the GPUs each rank sees
  (for example with ``ONEAPI_DEVICE_SELECTOR``) to assign them.
* An MPI library built with GPU support can create a context on the first GPU
  each rank sees while MPI initializes, before OpenSn selects a GPU. When
  ranks share GPUs, every rank then also holds memory on the node's first GPU.
  MPICH, for example, does so to register its host buffers. Restricting the
  GPUs each rank sees avoids this, and so does setting
  ``MPIR_CVAR_ENABLE_GPU_REGISTER=0`` for MPICH.

.. note::

   GPU acceleration is mainly a sweep-execution choice, not a different problem
   formulation. The same mesh, materials, groupsets, sources, and solver setup
   should usually be validated on CPU first.

GPU-aware MPI
-------------

With ``-DOPENSN_USE_GPU_AWARE_MPI=ON``, AAH sweeps on the GPU pass device
buffers directly to MPI instead of copying them through host memory. The MPI
library must support the GPU backend. CMake checks the library when it can:

* Cray MPICH: CMake links the GPU Transport Layer library (``libmpi_gtl_hsa``
  for HIP, ``libmpi_gtl_cuda`` for CUDA). Set ``MPICH_GPU_SUPPORT_ENABLED=1``
  when running. Without it, Cray MPICH treats device pointers as host memory.
* MPICH: configuration fails if the library was configured without the backend,
  for example ``--without-cuda``. If it was configured without
  ``--with-cuda``, ``--with-hip``, or ``--with-ze``, CMake reports that support
  could not be verified, since MPICH then enables the backend only if its
  configure step detected the GPU runtime.
* Open MPI: configuration fails if ``mpi-ext.h`` defines
  ``MPIX_CUDA_AWARE_SUPPORT`` (CUDA) or ``MPIX_ROCM_AWARE_SUPPORT`` (HIP) as 0,
  which means the library was built without that support. CMake warns if the
  macro is not defined, as in Open MPI versions without the extension.

On AMD MI300A nodes with Cray MPICH, GPU-aware AAH sweeps have shown
intermittent GPU memory access faults with Cray MPICH's intra-node GPU IPC
path. Set ``MPICH_GPU_IPC_ENABLED=0`` to use its other intra-node path. In
OpenSn scaling runs this cost no measurable time.

Typical execution pattern
=========================

Most inputs follow the same high-level structure:

1. Create or import the mesh.
2. Create materials and cross sections.
3. Define sources, boundaries, quadratures, and groupsets.
4. Construct the problem object.
5. Construct the solver object.
6. Call ``Initialize()``.
7. Run the solve with ``Execute()`` or an explicit transient time loop.
8. Create field functions or other outputs from the completed state.

For steady-state and eigenvalue problems, a typical driver looks like:

.. code-block:: python

   # This groupset covers all energy groups, numbered from zero.
   num_groups = groupset["groups_from_to"][1] + 1
   phys = opensn.DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=num_groups,
       groupsets=[groupset],
       xs_map=xs_map,
       boundary_conditions=boundary_conditions,
   )

   solver = opensn.SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

For transient problems, the driver can either call ``Execute()`` directly or
advance the state manually:

.. important::

   Transient problems require ``save_angular_flux=True`` in the problem
   options when constructing ``phys``. After creating ``phys``, call
   ``phys.SetTimeDependentMode()`` before constructing the
   :py:class:`pyopensn.solver.TransientSolver`.

.. code-block:: python

   phys.SetTimeDependentMode()
   dt = 1.0e-3
   stop_time = 1.0e-1
   num_steps = round(stop_time / dt)  # Fixed steps starting at time zero.

   solver = opensn.TransientSolver(
       problem=phys,
       dt=dt,
       stop_time=stop_time,
       initial_state="zero",
   )

   solver.Initialize()

   for _ in range(num_steps):
       solver.Advance()

Using ``Advance()`` is useful when the input needs to update sources,
boundaries, or other problem data between timesteps.

.. note::

   A manual transient loop is the right choice whenever the input itself needs
   to decide what happens next. If the timestep sequence is fixed and no
   mid-run updates are needed, ``Execute()`` is usually simpler.

Inputs, logs, and outputs
=========================

OpenSn writes its progress and solver information to standard output. During
development, this is often enough. For production runs, users commonly redirect
the output to a log file:

.. code-block:: bash

   mpirun -np 16 build/python/opensn -i problem.py > problem.log 2>&1

Field functions are created from the solved state after ``Execute()`` or after
an ``Advance()`` step in a transient loop. These can then be exported for
visualization or evaluated with post-processing utilities.

For transient problems, this includes the same requirement:
``save_angular_flux=True`` must already be enabled on the problem.

Restart and repeated solves
===========================

Because inputs are Python, it is natural to run multiple solves in one script.
For example, a study may:

1. Build a base problem.
2. Solve it once.
3. Change a source or boundary condition.
4. Solve again.

This is a useful pattern, but it is important to keep the distinction between
problem data and solver state in mind. If a change affects the transport
operator or the source model, make sure the updated objects are set on the
problem before executing the next solve.

Practical advice
================

- Start by running new inputs on one MPI rank before scaling out.
- Keep the first version of an input simple: one groupset, one source, one
  field-function export.
- Add post-processing only after the physics setup and solver convergence are
  behaving as expected.
- For transient problems, call ``Update()`` on existing field functions after
  each completed ``Advance()`` step, or create fresh field functions from the
  current state. Do not assume an earlier field-function object updates
  automatically.

.. note::

   Many input problems that look like "solver failures" are actually setup
   problems: wrong material assignment, an inconsistent groupset split, a
   missing source, or a boundary condition that does not match the intended
   physical model. When in doubt, reduce the input to the smallest case that
   should still work and build back up from there.
