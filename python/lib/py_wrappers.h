// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

#include <pybind11/pybind11.h>

#include "framework/parameters/parameter_block.h"

namespace py = pybind11;

namespace opensn
{

// Converter
// ---------

/// Convert a C++ vector to a Python memoryview.
template <typename T>
py::memoryview convert_vector(const std::vector<T>& vec)
{
  return py::memoryview::from_buffer(const_cast<T*>(vec.data()), {vec.size()}, {sizeof(T)}, true);
}

/// Translate a Python dictionary into a ParameterBlock.
ParameterBlock kwargs_to_param_block(const py::kwargs& params);

// Module wrappers
// ---------------

/// Wrap the angular quadrature components of OpenSn.
void py_aquad(py::module& pyopensn);
void WrapQuadrature(py::module& aquad);
void WrapProductQuadrature(py::module& aquad);
void WrapCurvilinearQuadrature(py::module& aquad);
void WrapSLDFESQuadrature(py::module& aquad);

/// Wrap the field function components of OpenSn (unfinished).
void py_ffunc(py::module& pyopensn);
void WrapFieldFunction(py::module& ffunc);
void WrapFieldFunctionGridBased(py::module& ffunc);
void WrapFieldFunctionInterpolation(py::module& ffunc);

/// Wrap the logical volume components of OpenSn (unfinished).
void py_logvol(py::module& pyopensn);
void WrapLogicalVolume(py::module& logvol);

// Wrap the math components of OpenSn
void py_math(py::module& pyopensn);
void WrapYlm(py::module& math);
void WrapVector3(py::module& math);
void WrapFunctors(py::module& math);

/// Wrap the mesh components of OpenSn.
void py_mesh(py::module& pyopensn);
void WrapMesh(py::module& mesh);
void WrapMeshGenerator(py::module& mesh);
void WrapGraphPartitioner(py::module& mesh);

/// Wrap the solver components of OpenSn (unfinished).
void py_solver(py::module& pyopensn);
void wrap_solver(py::module& slv);

/// Wrap the source components of OpenSn (unfinished).
void py_source(py::module& pyopensn);
void WrapPointSource(py::module& src);
void WrapVolumetricSource(py::module& src);

/// Wrap the cross section components of OpenSn.
void py_xs(py::module& pyopensn);
void WrapMultiGroupXS(py::module& xs);

} // namespace opensn
