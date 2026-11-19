[![Regression Tests](https://github.com/Open-Sn/opensn/actions/workflows/regression.yaml/badge.svg)](https://github.com/Open-Sn/opensn/actions/workflows/regression.yaml)
[![Weekly Tests](https://github.com/Open-Sn/opensn/actions/workflows/weekly.yaml/badge.svg)](https://github.com/Open-Sn/opensn/actions/workflows/weekly.yaml)
[![Nightly Builds](https://github.com/Open-Sn/opensn/actions/workflows/builds.yaml/badge.svg)](https://github.com/Open-Sn/opensn/actions/workflows/builds.yaml)
[![Coverage Status](https://coveralls.io/repos/github/Open-Sn/opensn/badge.svg?branch=main)](https://coveralls.io/github/Open-Sn/opensn?branch=main)
[![License: MIT](https://img.shields.io/badge/license-MIT-008b8b)](https://github.com/Open-Sn/opensn/blob/main/LICENSE)
[![Documentation](https://img.shields.io/github/actions/workflow/status/Open-Sn/opensn/gh-pages.yml?branch=main&label=docs)](https://open-sn.github.io/opensn/)

<img src="doc/assets/hero.png" alt="OpenSn — Radiation transport across every scale. C5G7 Thermal Flux and Tokamak Neutron Flux." width="100%">

**Open-source, massively parallel radiation transport.** OpenSn solves the
multigroup linear Boltzmann transport equation with the discrete-ordinates
($S_N$) method on structured and unstructured meshes. A modern C++ engine and
Python interface bring reactor physics, shielding, and charged-particle
transport into one simulation workflow—from a laptop to a supercomputer.

[Documentation](https://open-sn.github.io/opensn/) ·
[Installation](https://open-sn.github.io/opensn/install/install.html) ·
[Tutorials](https://open-sn.github.io/opensn/tutorials/index.html) ·
[Examples](https://github.com/Open-Sn/OpenSnExamples) ·
[Discussions](https://github.com/Open-Sn/opensn/discussions)

## See what OpenSn can do...

Below are some examples of the types of problems that OpenSn can solve.

<div>
<a href="https://github.com/Open-Sn/OpenSnExamples/tree/main/CASL_VERA_2D"><img src="doc/assets/gallery-c5g7.png" alt="Explore Reactor Physics: relative thermal-group scalar flux in a multigroup reactor criticality calculation." width="50%" align="top"></a><a href="https://github.com/Open-Sn/OpenSnExamples/tree/main/Coiled_Wire_Electron_Photon"><img src="doc/assets/gallery-wire.png" alt="Explore Charged Particle Transport: OpenSn CSDA energy deposition in a coiled wire." width="50%" align="top"></a>
</div>
<div>
<a href="https://github.com/Open-Sn/OpenSnExamples/tree/main/Urban_Source"><img src="doc/assets/gallery-urban.png" alt="Explore Shielding and Streaming: OpenSn 25-group effective dose rate from a Cs-137 source in the urban-city example." width="50%" align="top"></a><a href="https://github.com/Open-Sn/OpenSnExamples/tree/main/Kobayashi_Time_Dependent"><img src="doc/assets/gallery-transient.gif" alt="Explore Time Dependent Transport: animated Kobayashi dog-leg scalar flux from 2 to 200 seconds." width="50%" align="top"></a>
</div>
<div>
<a href="https://github.com/Open-Sn/OpenSnExamples/tree/main/OpenFOAM_Electron_Converter"><img src="doc/assets/gallery-multiphysics.png" alt="Explore Coupled Multiphysics: OpenSn energy deposition in a cutaway water-cooled converter coupled to OpenFOAM." width="50%" align="top"></a><a href="https://github.com/Open-Sn/OpenSnExamples/tree/main/Tokamak_Fusion_Engineering"><img src="doc/assets/gallery-fusion.png" alt="Explore Fusion Engineering: total neutron flux in the OpenSn tokamak shielding problem." width="50%" align="top"></a>
</div>

## Capabilities

- **Particle transport:** multigroup neutron and photon transport, coupled
  particle groups, and steady-state CSDA electron/positron transport.
- **Problem types:** steady-state and transient fixed-source problems,
  criticality ($k$-eigenvalue) calculations, and steady-state adjoint responses
  and cross-section sensitivities.
- **Geometry:** 1D, 2D, and 3D discontinuous finite elements on polygonal and
  polyhedral meshes; spherical and axisymmetric curvilinear transport.
- **Angular resolution:** product quadratures and locally refined spherical
  finite-element quadratures; uncollided/first-collision treatments for
  supported point-source problems.
- **Solver acceleration:** Krylov iterative algorithms, diffusion synthetic
  acceleration, SMM, and CMFD.
- **Parallel execution:** MPI domain decomposition and GPU sweep backends for
  CUDA, HIP, and SYCL.
- **Modeling and analysis:** Python scripting and C++ APIs, native/OpenMC/CEPXS
  cross-section readers, restart files, integral responses, and VTK output.

See the [user guide](https://open-sn.github.io/opensn/userguide/index.html) for details.

## Get started

1. **Install OpenSn** using the [installation guide](https://open-sn.github.io/opensn/install/install.html).
2. **Run a small problem** from the [tutorials](https://open-sn.github.io/opensn/tutorials/index.html).
3. **Build on a benchmark** from [OpenSnExamples](https://github.com/Open-Sn/OpenSnExamples),
   including CASL VERA, shielding, and source-driven problems.

OpenSn is developed at Texas A&M University with collaborators and contributors
from the wider radiation-transport community. Questions and modeling discussions
are welcome in [GitHub Discussions](https://github.com/Open-Sn/opensn/discussions).
For bugs, use [Issues](https://github.com/Open-Sn/opensn/issues); for contributions,
start with the [developer workflow](https://github.com/Open-Sn/opensn/blob/main/doc/source/devguide/workflow.rst).

## Cite OpenSn

If OpenSn supports your research, please cite
[Andrs et al., *OpenSn: A Massively Parallel, Open-Source Simulation Environment for Discrete Ordinates Radiation Transport*](https://doi.org/10.1016/j.cpc.2025.110001),
*Computer Physics Communications* **320**, 110001 (2026).
OpenSn is distributed under the [MIT license](https://github.com/Open-Sn/opensn/blob/main/LICENSE).
