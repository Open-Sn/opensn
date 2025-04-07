// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/parameter_block.h"
#include "framework/math/vector.h"
#include <pybind11/pybind11.h>
#include <vector>

namespace py = pybind11;

namespace opensn
{

// Converter

/// Convert a C++ vector to a Python memoryview.
template <typename T>
py::memoryview
convert_vector(const std::vector<T>& vec)
{
  return py::memoryview::from_buffer(const_cast<T*>(vec.data()), {vec.size()}, {sizeof(T)}, true);
}

/// Convert an OpenSn vector to a Python memoryview.
template <typename T>
py::memoryview
convert_vector(const Vector<T>& vec)
{
  return py::memoryview::from_buffer(const_cast<T*>(vec.data()), {vec.size()}, {sizeof(T)}, true);
}

/// Translate a Python dictionary into a ParameterBlock.
ParameterBlock kwargs_to_param_block(const py::kwargs& params);

// Module wrappers

/// Wrap the angular quadrature components of OpenSn.
void py_aquad(py::module& pyopensn);
void WrapQuadraturePointPhiTheta(py::module& aquad);
void WrapQuadrature(py::module& aquad);
void WrapProductQuadrature(py::module& aquad);
void WrapCurvilinearQuadrature(py::module& aquad);
void WrapSLDFESQuadrature(py::module& aquad);

// Wrap the diffusion components of OpenSn
void py_diffusion(py::module& pyopensn);
void WrapDiffusion(py::module& diffusion);

/// Wrap the field function components of OpenSn.
void py_ffunc(py::module& pyopensn);
void WrapFieldFunction(py::module& ffunc);
void WrapFieldFunctionGridBased(py::module& ffunc);
void WrapFieldFunctionInterpolation(py::module& ffunc);

/// Wrap the logical volume components of OpenSn.
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

/// Wrap the response components of OpenSn.
void py_response(py::module& pyopensn);
void WrapResEval(py::module& response);

/// Wrap the post-processing components of OpenSn (unfinished).
void py_post(py::module& pyopensn);
void WrapPostProcessor(py::module& post);
void WrapPrinter(py::module& post);

/// Wrap the settings components of OpenSn
void py_settings(py::module& pyopensn);

/// Wrap the solver components of OpenSn (unfinshed).
void py_solver(py::module& pyopensn);
void WrapSolver(py::module& slv);
void WrapLBS(py::module& slv);
void WrapSteadyState(py::module& slv);
void WrapNLKEigen(py::module& slv);
void WrapPIteration(py::module& slv);
void WrapPRK(py::module& slv);

// Wrap the diffusion components of OpenSn
void py_diffusion(py::module& pyopensn);
void WrapDiffusion(py::module& diffusion);

/// Wrap the source components of OpenSn.
void py_source(py::module& pyopensn);
void WrapPointSource(py::module& src);
void WrapVolumetricSource(py::module& src);

/// Wrap the cross section components of OpenSn.
void py_xs(py::module& pyopensn);
void WrapMultiGroupXS(py::module& xs);

} // namespace opensn
