# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

# The oldest version of each dependency that OpenSn supports. CMakeLists.txt, the installed
# OpenSnConfig.cmake, tools/dependencies, and the documentation read these values, and
# tools/developer/check-spack-floors.cmake keeps the Spack recipe's floors at or above them.
# CMake does not check the HDF5, VTK, Caliper, and Python floors yet; only the documentation and
# the Spack check use them. doc/source/conf.py parses this file, so keep each floor a one-line
# set().

set(OPENSN_MIN_MPICPP_LITE_VERSION 2.7.1)
set(OPENSN_MIN_BOOST_VERSION 1.86)
set(OPENSN_MIN_PETSC_VERSION 3.17)
set(OPENSN_MIN_HDF5_VERSION 1.14)
set(OPENSN_MIN_VTK_VERSION 9.3)
set(OPENSN_MIN_CALIPER_VERSION 2.11)
set(OPENSN_MIN_PRECICE_VERSION 3.4)
set(OPENSN_MIN_CUDA_VERSION 12.0)
set(OPENSN_MIN_PYTHON_VERSION 3.9)
