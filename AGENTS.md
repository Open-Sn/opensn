# OpenSn Agent Guide

This file applies to the entire repository. Use it with the developer and
installation documentation under `doc/source/`. A more specific `AGENTS.md`
takes precedence within its subtree.

## Priorities

- Treat OpenSn as production scientific software. Numerical correctness,
  parallel correctness, and reproducibility are release requirements.
- Verify reports and review comments against current code. Trace callers,
  ownership, MPI participation, and CPU/GPU variants before editing.
- Make the smallest coherent change. Avoid unrelated cleanup, broad renaming,
  generated output, and formatting untouched files.
- Inspect `git status` and the relevant diff first. Preserve user changes and
  untracked files. Do not commit, push, rebase, or change branches unless asked.
- Use repository scripts and existing test infrastructure. Report exact
  validation and anything that could not be run.

### Numerical correctness

Implementation, documentation, and tests must agree on equations,
discretizations, units, signs, normalization, indexing, boundary conditions,
residuals, and convergence criteria.

Compilation, convergence, or agreement between two paths does not establish
correctness. Anchor affected results in an analytic or manufactured solution,
conservation law, mathematical limit, published benchmark, or independently
generated reference. Quantify error with tolerances justified by
discretization, iteration, and floating-point effects.

Do not conceal a discrepancy by loosening a tolerance, replacing a gold value
with output from the changed implementation, suppressing a residual, or
reducing coverage. Establish which result is correct. Performance changes must
preserve the mathematical solution within justified tolerances. Make any
accuracy-changing approximation explicit, documented, optional when
appropriate, and independently validated.

For solver changes, check conservation and relevant invariants as well as the
primary answer. Distinguish convergence of an algebraic iteration from accuracy
with respect to the underlying problem. If correctness cannot be established,
state the uncertainty and the validation still required.

## Repository and environment

- `framework/`: mesh, math, MPI, runtime, logging, and shared infrastructure.
- `modules/`: physics solvers, including linear Boltzmann transport.
- `python/`: console application and Python bindings.
- `pyopensn/`: Python module package.
- `test/unit/`: GoogleTest unit tests.
- `test/python/`: regression inputs and `tests.json` specifications.
- `doc/source/`: Sphinx manuals and executable tutorials.
- `tools/developer/`: sanitizer and clang-tidy tooling.
- `external/`: vendored code; edit only when the task concerns that dependency.

Find current compiler, dependency, and tool requirements in
`doc/source/install/`, the root build configuration, and `tools/dependencies/`.
Do not copy version requirements into this file.

Use one consistent compiler/MPI/Python/dependency stack. The CMake compiler,
MPI wrappers, PETSc, HDF5, Python extension, and `mpi4py` must be ABI compatible.
Use separate build directories when changing compiler, MPI, build type, or GPU
backend. Diagnose stale caches and mixed environments before changing source.

For reproducible failures and measurements, record the commit, compiler, MPI,
Python environment, CMake options, build type, hardware/backend, ranks, threads,
mesh, cross sections, quadrature, groups, solver settings, tolerances, and seeds.
Do not commit machine paths, build trees, profiling databases, credentials, or
site-specific launcher settings.

## Build and package

Normal developer build with console and Python module:

```sh
cmake -S . -B build -DOPENSN_WITH_PYTHON_MODULE=ON
cmake --build build -j
```

Debug build:

```sh
cmake -S . -B build-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DOPENSN_WITH_PYTHON_MODULE=ON
cmake --build build-debug -j
```

Sanitizer build and representative regressions:

```sh
cmake -DOPENSN_WITH_PYTHON_MODULE=ON --preset clang+debug+sanitizer
cmake --build --preset clang+debug+sanitizer -j
export LSAN_OPTIONS="suppressions=$PWD/tools/developer/lsan.supp"
test/run_tests -d test/python -j 4 -v 1 -w 3 \
  --exe build-debug-sanitizer/python/opensn
```

Build the affected target while iterating, then complete the applicable build.
Test a clean configure when CMake, dependencies, feature options, install rules,
or source registration changes.

### GPU builds

Only one GPU backend may be enabled per build. For CUDA, use the architecture
appropriate to the machine. Set `CUDA_ARCH` from the active toolchain or site
configuration before running:

```sh
cmake -S . -B build-cuda \
  -DOPENSN_WITH_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH" \
  -DOPENSN_WITH_PYTHON_MODULE=ON
cmake --build build-cuda -j
```

HIP and SYCL use their corresponding `OPENSN_WITH_*` options and architecture
settings. Never guess available hardware. CPU-only configurations must remain
free of GPU requirements, and optional Python/module builds must remain valid.

### Build-system and dependency changes

- Keep includes, definitions, options, and libraries at the narrowest CMake
  target scope. Avoid global flags and leaked build-tree paths.
- Confirm new or moved files enter the intended target and compilation database.
- Update install/export rules for public headers, libraries, and targets.
- Justify new dependencies, preserve supported discovery workflows, verify
  license compatibility, and consider ABI and static/shared variants.
- Keep feature detection explicit and configuration failures actionable.
- Do not weaken repository-wide warnings to accommodate one target or toolchain.

For public packaging changes, install to a fresh prefix and exercise it outside
the source and build trees:

```sh
cmake --install build --prefix /tmp/opensn-install
```

Test C++ exports with a minimal external `find_package(OpenSn CONFIG REQUIRED)`
consumer. Test Python packaging in a fresh environment without the source-tree
`PYTHONPATH`. Verify runtime discovery of OpenSn and dependent libraries; an
in-place build alone is insufficient.

## Run, debug, and analyze

Prefer the regression runner because `tests.json` supplies the working
directory, ranks, arguments, and checks:

```sh
test/run_tests --exe build/python/opensn \
  -t test/python/path/to/test.py -j 4 -v 1
```

Direct console run:

```sh
build/python/opensn -i path/to/input.py
```

The console and `pyopensn` module are separate modes; do not import `pyopensn`
inside the console. Exercise module-visible changes in module mode:

```sh
export PYTHONPATH="$PWD/build${PYTHONPATH:+:$PYTHONPATH}"
test/run_tests --engine=module -t path/to/module_test.py -j 4 -v 1
```

Use a Debug build with `gdb` or `lldb` for control flow and state. Use the
sanitizer preset for invalid accesses, lifetime errors, leaks, and undefined
behavior. Also reproduce optimization-sensitive failures with `RelWithDebInfo`.
Suppress a sanitizer or clang-tidy finding only after demonstrating a tool or
third-party false positive, and scope the suppression narrowly.

### Parallel and accelerator analysis

- All ranks must enter collectives in the same order. Check root/nonroot buffer
  sizes, count types, tags, requests, communicator lifetime, and zero-work ranks.
- Reproduce MPI defects at the smallest rank count that exercises them, then at
  the rank count required by the regression.
- For hangs, inspect unmatched communication, divergent collectives, unresolved
  sweep dependencies, and incomplete GPU events.
- Test shared-state changes with `OPENSN_NUM_THREADS=1` and a representative
  multithreaded value. Check races, nested parallelism, thread-local lifetime,
  binding, and oversubscription.
- For GPU code, inspect host/device ownership, transfer sizes, launch bounds,
  partial work groups, empty buffers, asynchronous lifetime, synchronization,
  error propagation, and CPU/GPU numerical parity.
- For `mpi4py`, use the MPI implementation that built OpenSn and `mpi4py`; check
  communicator ownership and initialization/finalization ordering.

### Memory, scalability, and resources

- Derive storage and communication cost in cells, faces, groups, angles, and MPI
  ranks. Avoid replicated `O(P^2)` data unless required and demonstrably bounded.
- Separate setup-only data from execution and rebuild state. Define repeated-call
  behavior and keep reads valid for the documented lifetime.
- Pair owned MPI, PETSc, HDF5, and GPU handles with release after last use,
  respecting initialization/finalization boundaries. Prefer RAII.
- Validate counts, dimensions, IDs, offsets, and buffer-size arithmetic before
  allocation or indexing. Use fixed-width persisted types and checked
  conversions; reject malformed or inconsistent input clearly.
- Avoid references, pointers, spans, or device views that outlive containers or
  asynchronous operations. Do not depend on native width, byte or container
  order, undefined evaluation, or uninitialized data.
- Check empty partitions, minimum sizes, cycles, reflecting boundaries, multiple
- Build validated temporary state before publishing it so failures leave owned
  resources releasable.

Profile before optimizing. Separate setup, graph construction, sweeps, solver
work, reductions, transfers, I/O, and finalization. Record hardware placement,
build options, problem dimensions, convergence, samples, variation, component
times, peak memory, and communication. Compare equivalent converged solutions;
do not infer scalability from one small run or use fragile wall-clock gates in
routine correctness tests.

### Compatibility, files, and diagnostics

Treat documented APIs, defaults, enum/string values, output keys, restart data,
meshes, cross sections, HDF5 data, and exported results as public contracts.
Search all consumers before changing them. For intentional breaks, update
callers and documentation together and provide a deprecation path when
reasonable. Never silently change a default that affects physics or convergence.

For persistent formats, document units, shapes, ordering, numeric types,
normalization, ownership, and versioning. Validate metadata before allocation.
Test round trips, empty data, incomplete or incompatible files, and multiple MPI
decompositions. Avoid native-width or container-order assumptions and provide a
migration path when compatibility cannot be preserved.

Errors should identify the operation and relevant rank, object, index, file, or
parameter without dumping datasets. Keep global summaries on one rank and
rank-local details rank-labeled. Avoid logging in hot loops and never change
collective participation merely to produce diagnostics.

## Tests

Run the narrowest relevant test while developing, then broaden validation in
proportion to the change.

| Test | Use |
| --- | --- |
| Unit | Algorithms, data structures, validation, ownership, and lifecycle |
| Regression | Complete OpenSn inputs and numerical behavior |
| MPI | Collectives, partitioning, communication, and distributed state |
| GPU | Backend behavior and CPU/GPU equivalence |
| Sanitizer | Memory, lifetime, bounds, leaks, and undefined behavior |
| Restart | Persistence and reconstruction of solver state |
| Tutorial | Executable end-user workflows |
| Benchmark | Runtime, memory, communication, and scaling outside routine correctness gates |

### Running tests

```sh
# Unit tests
build/test/opensn-unit
mpirun -np 2 build/test/opensn-unit
mpirun -np 4 build/test/opensn-unit
build/test/opensn-unit --gtest_filter='SuiteName.TestName'

# Selected regression
test/run_tests --exe build/python/opensn -t path/to/test.py -j 4 -v 1

# Normal CPU regression set
test/run_tests -d test/python -j 8 -v 1 -w 3 \
  --exe build/python/opensn

# GPU regressions, with a configured GPU build and available hardware
test/run_tests --gpu -d test/python -j 8 -v 1 \
  --exe build-cuda/python/opensn
```

`-j` is the total CPU-slot budget; MPI tests consume their requested ranks. Do
not oversubscribe shared resources.

### Creating unit tests

Add C++ tests under the matching `test/unit/` subtree and follow nearby
GoogleTest organization. Test one public contract with small deterministic
fixtures. Cover nominal, empty/minimum, boundary, invalid, repeated-call,
mutation, and failure-cleanup cases as applicable. Test copy/move behavior only
when the type supports it.

Use exact comparisons for discrete results and derived tolerances for floating
point. Prefer `ASSERT_*` only when continuing is unsafe; otherwise use
`EXPECT_*` for better diagnostics. Verify documented exception types without
depending on incidental wording. Do not test private implementation details,
duplicate the algorithm in the test, or depend on test order, timing, ambient
files, or leaked global state.

A bug fix should have a focused test that fails for the defect when practical.
If the defect requires MPI, GPU hardware, or complete solver setup, use the
smallest appropriate regression instead.

### Creating regression tests

Place inputs beside related tests under `test/python/` and register them in the
nearest `tests.json`. Set appropriate `num_procs`, `weight_class`, checks,
tolerances, and a unique `outfileprefix` when reusing an input.

State the mathematical claim and reference in the input. Prefer, in order:

1. exact or semi-analytic solutions;
2. conservation identities;
3. manufactured solutions with independently evaluated sources;
4. limits or invariants such as zero source, symmetry, positivity, or infinite
   medium behavior;
5. published benchmarks with citations; or
6. independently generated higher-fidelity results with reproducible settings.

A value generated only by the code under test is a baseline, not proof of
correctness. Equivalence tests for acceleration, sweeps, MPI, rebuilds, and GPU
paths should share an independently anchored reference because common defects
can agree.

Use a compact but realistic model that activates the production path, including
relevant heterogeneity, boundaries, energy coupling, mesh type, and partitioning.
Do not use a trivial model that bypasses the feature, or make routine CI tests
large when scale is not essential.

Derive acceptance tolerances from spatial, angular, energy, temporal, iterative,
and floating-point error. Keep solver tolerances tighter than the acceptance
tolerance. Use absolute and relative tolerance near zero and allow normal
reduction/backend variation without permitting the original defect. Demonstrate
expected convergence across resolutions when that is the mathematical claim.

Print stable, unique metrics from one rank. Prefer numeric `KeyValuePair` or
`FloatCompare` checks for physical quantities, errors, norms, balances, rates,
or selected points. `ErrorCode` alone rarely establishes scientific correctness.
Avoid large output dumps and clean temporary files and global state.

Before accepting a test, confirm it fails for the defect or a deliberate
perturbation, verify the reference independently, repeat it for determinism,
check meaningful serial/MPI/CPU/GPU variants, and ensure missing output or early
exit cannot pass. Re-derive any changed gold value and document its provenance.

For public APIs, test defaults, types, ownership, invalid inputs, and collective
requirements. For restart, compare uninterrupted and resumed results and reject
bad data. Compare MPI decompositions within reduction tolerance, including
zero-work ranks when supported. Compare GPU paths with an anchored CPU result.
Keep performance benchmarks separate from correctness gates.

## Style and static analysis

Follow `doc/source/devguide/coding_standard.rst`,
`doc/source/devguide/workflow.rst`, and the root `.clang-format`, `.clang-tidy`,
and `.flake8`. Match nearby naming, ownership, exception, and header conventions.
New C++ source and headers require the repository SPDX header.

Format changed tracked C/C++ and GPU files, and run `clang-format -i` directly
on new untracked source files:

```sh
git diff --name-only -z --diff-filter=ACMR HEAD -- \
  '*.c' '*.cc' '*.cpp' '*.cxx' '*.h' '*.hh' '*.hpp' '*.cu' \
  ':(exclude)doc/**' ':(exclude)resources/**' \
  ':(exclude)tutorials/**' ':(exclude)external/**' \
  | xargs -0 -r clang-format -i
```

Check formatting by replacing `-i` with `--dry-run --Werror`.

clang-tidy requires `build/compile_commands.json` and treats findings as errors:

```sh
tools/developer/run-clang-tidy.sh path/to/changed_file.cc
tools/developer/run-clang-tidy.sh  # repository-wide when warranted
```

For CUDA, use the repository wrapper with the configured build, toolkit, and
machine architecture. Set `CUDA_ROOT` and `CUDA_CLANG_ARCH` from the active
environment:

```sh
tools/developer/run-cuda-clang-tidy.sh \
  --build-dir build-cuda \
  --cuda-path "$CUDA_ROOT" \
  --cuda-arch "$CUDA_CLANG_ARCH" \
  path/to/changed_file.cu
```

Python CI uses:

```sh
flake8 . --count --show-source --statistics
```

Apply checks to changed files first. Do not rewrite unrelated files to resolve
pre-existing repository-wide findings.

## Documentation

Documentation is part of a feature. Update it when behavior, interfaces,
defaults, ranges, units, lifecycle, collectivity, normalization, output, or
supported cases change. Document current behavior, not plans.

Binding docstrings in `python/lib/` and `doc/source/pyapi/` define the Python API;
public-header Doxygen and `doc/source/capi/` cover C++. User, theory, tutorial,
and developer material belongs in the corresponding `doc/source/` subtree.
Significant solver features often need API documentation, user guidance, theory,
and a tested example.

### API documentation

Follow `doc/source/devguide/py_bindings.rst` for NumPy-style Python docstrings
and `doc/source/devguide/doxygen.rst` for C++. Document purpose, parameters,
types, units, defaults, valid ranges, alternatives, return ownership/shape,
side effects, lifecycle, exceptions, collective MPI requirements, and backend
restrictions. Signatures, `py::arg` defaults, validation, and prose must agree.

Add public Python objects to the appropriate autosummary list in
`doc/source/pyapi/index.rst`. Do not edit generated autosummary pages. Put
Doxygen contracts on header declarations and describe ownership, lifetime,
valid states, complexity, and thread/MPI requirements where relevant.

### User guide, theory, and tutorials

User-guide material should explain when to use a feature, supported cases,
complete setup, option interactions, units, accuracy/performance implications,
execution order, common mistakes, and related API/theory/tutorial pages. Keep
examples executable and include every required argument.

Theory must define symbols, measures, indices, units, signs, norms,
normalization, and assumptions. Distinguish general theory from OpenSn's
implemented discrete convention. Trace important quantities to the algebraic
solve and document residuals, stopping tests, updates, safeguards, and skipped
operations. Check dimensions, group directions, moments, adjoint transposes,
interface signs, limits, and code agreement. Cite nontrivial methods and explain
implemented approximations.

Tutorials should teach one complete task with a small realistic model, explain
the expected physical result, and include useful postprocessing. Anchor results
in mathematics, conservation, a benchmark, or clearly explained qualitative
behavior. Register executable tutorials and avoid unnecessary generated assets.

Use consistent terminology and Sphinx cross-reference roles. Add pages to the
appropriate toctree. Prefer `literalinclude` for snippets that must remain
synchronized with tested files. Give figures labels, units, captions, and
accessible descriptions.

### Documentation validation

```sh
python -m pip install -r doc/requirements.txt
export PYTHONPATH="$PWD/build${PYTHONPATH:+:$PYTHONPATH}"
make -C doc html

test/run_tests -d doc/source/tutorials -j 4 -v 1 -w 3 \
  --exe build/python/opensn --engine=jupyter
```

Use `SPHINXOPTS='-W --keep-going'` for a strict final build when the tree permits
it. Run/import changed examples: Sphinx cannot prove code executes or equations
match implementation. Verify signatures, defaults, units, MPI collectivity,
supported cases, theory conventions, references, and related tests.

## Licensing, generated assets, and CI

Use the SPDX form found in neighboring source files. Preserve third-party
notices and verify redistribution rights for code, equations, figures, meshes,
cross sections, and benchmark data. Record scientific asset source, units,
normalization, transformations, license, and generation procedure. Prefer small
reviewable fixtures and reproducible generators over opaque binaries.

Edit a generator rather than generated API pages, Doxygen output, meshes, or
reference results. If a generated artifact is tracked, verify deterministic
regeneration and expected-only changes.

Keep CI changes scoped and equivalent to developer commands. Preserve compiler,
backend, rank, sanitizer, logging, and artifact coverage. Make skips explicit,
avoid shared mutable output, validate quoting/configuration, and do not turn a
failure into success to make a patch green.

## Review and contribution

Review the actual diff and call paths; do not assume a report or historical
implementation is correct. Inspect lifecycle and ownership, failure cleanup,
edge cases, MPI/thread/GPU paths, numerical definitions and conservation,
production-scale complexity, public compatibility, and test sensitivity.

Report findings only with a concrete trigger and current-code evidence. Give the
location, preconditions, impact, and smallest reasonable fix. Distinguish
confirmed defects from platform- or scale-dependent risks. After fixing, rerun
analysis and relevant tests; for scalability findings, state old and new
asymptotic costs and expected hot-path effects.

Development belongs on a feature branch. Follow
`doc/source/devguide/workflow.rst`: rebase onto `main`; do not merge `main` into
the feature branch. Keep commits focused with short imperative subjects.

Before handing work back:

1. Inspect `git diff`, `git diff --check`, and status for unrelated files.
2. Verify SPDX, licensing, provenance, and generated-file handling.
3. Run formatting, clang-tidy, and Python linting as applicable.
4. Build the affected configuration; clean-configure build-system changes.
5. Run focused tests plus relevant MPI, threading, sanitizer, GPU, restart,
   binding, tutorial, and documentation validation.
6. Check API/file compatibility and installed-tree use when affected.
7. Independently validate equations, references, tolerances, and gold values.
8. Measure representative performance for hot-path or complexity changes.
9. Report behavior, rationale, exact validation, and remaining platform, scale,
   backend, or compatibility limits.
