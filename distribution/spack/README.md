# OpenSn Spack Repository

This directory contains the local Spack package repository for OpenSn. It
provides recipes for:

- `opensn`
- `mpicpp-lite`
- a PETSc override used to support the PETSc/Hypre constraints required by
  OpenSn

## Add the Repository

From an OpenSn checkout:

```bash
spack repo add $PWD/distribution/spack
```

Confirm that Spack sees the local package:

```bash
spack info opensn
```

## Install

Default CPU build:

```bash
spack install opensn@1.0.1
```

Build the current `main` branch:

```bash
spack install opensn@develop
```

Build with OpenSn's native optimized CMake configuration:

```bash
spack install opensn@1.0.1+native
```

Build the `pyopensn` module in addition to the Python-enabled console:

```bash
spack install opensn@1.0.1+python_module
```

Build with a GPU backend:

```bash
spack install opensn@1.0.1+cuda
spack install opensn@1.0.1+hip
spack install opensn@1.0.1+sycl sycl_arch=spir64_gen
```

CUDA, HIP, and SYCL are mutually exclusive.

## Supported Variants

OpenSn-specific variants:

- `+python` / `~python`: build the Python-enabled `opensn` console. Default:
  `+python`.
- `+python_module` / `~python_module`: build and install `pyopensn`. Default:
  `~python_module`. Requires `+python`.
- `+native` / `~native`: use OpenSn's `Native` CMake build type. Default:
  `~native`.
- `+cuda`, `+hip`, `+sycl`: enable one GPU backend. Default: all disabled.
- `sycl_arch=<value>`: SYCL target architecture. Default: `spir64_gen` (see [ICPX](https://github.com/intel/llvm/blob/sycl/sycl/doc/UsersManual.md) or [ACPP](https://github.com/AdaptiveCpp/AdaptiveCpp/blob/develop/doc/using-acpp.md)).
- `sycl_target_backend=<value>`: optional SYCL target backend value.
- `sycl_register_alloc_mode=<value>`: optional SYCL register allocation mode.

The package also inherits standard `CMakePackage` variants such as
`build_type`, `generator`, and `ipo`.

## Tests

Run unit and regression tests during installation:

```bash
spack install --test=root opensn@1.0.1
```

Run post-install tests:

```bash
spack test run opensn
spack test results
```

The post-install test support caches `test/run_tests`, `test/src`,
`test/python`, and `test/assets`, and installs a private copy of the unit-test
binary under `libexec/opensn/opensn-unit`.

## Locate the Executable

```bash
spack location -i opensn
ls "$(spack location -i opensn)/bin/opensn"
```

If multiple OpenSn installations exist, use the package hash:

```bash
spack find -l opensn
spack location -i /<hash>
```

## Development Workflow

Use a Spack environment when you want Spack to manage OpenSn's dependency
stack but keep the OpenSn source tree editable.

From the OpenSn checkout:

```bash
spack env create opensn-dev
spack env activate opensn-dev
spack repo add $PWD/distribution/spack
spack add opensn@develop
spack develop --path $PWD opensn@develop
spack concretize
spack install
```

After installation, load the environment and use the installed OpenSn
executable:

```bash
spack env activate opensn-dev
spack load opensn
opensn --version
```

For an edit-build-test loop, rebuild the developed package after changing the
source tree:

```bash
spack install opensn
```

If you want to configure and build manually while still using Spack-built
dependencies, activate the environment, install the dependencies, and point
CMake at Spack's view:

```bash
spack env activate opensn-dev
spack install --only dependencies opensn
spack env view enable
cmake -S . -B build -DCMAKE_PREFIX_PATH=$SPACK_ENV/.spack-env/view
cmake --build build
```

## Notes

- The OpenSn package models the dependencies required by the top-level
  `CMakeLists.txt`: MPI, `mpicpp-lite`, Boost, PETSc, HDF5, VTK, Caliper,
  GoogleTest, Python, and pybind11.
- PETSc is constrained to the OpenSn-required configuration: MPI-enabled,
  double precision, real scalars, 64-bit indices, HDF5, METIS, PTScotch, and
  Hypre.
- The local PETSc override is intentionally close to Spack's builtin PETSc
  recipe. Spack patch files are referenced by pinned URLs instead of copied
  into this repository.
- `+native` maps to `-DCMAKE_BUILD_TYPE=Native`, which enables OpenSn's native
  optimization flags when supported by the compiler.
