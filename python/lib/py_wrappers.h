// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/parameter_block.h"
#include <pybind11/pybind11.h>
#include <vector>

namespace py = pybind11;

namespace opensn
{

/// Convert a C++ vector to a Python memoryview.
template <typename T>
py::memoryview
convert_vector(const std::vector<T>& vec)
{
  return py::memoryview::from_buffer(const_cast<T*>(vec.data()), {vec.size()}, {sizeof(T)}, true);
}

/// Translate a Python dictionary into a ParameterBlock.
ParameterBlock kwargs_to_param_block(const py::kwargs& params);

/// Wrap the angular quadrature components of OpenSn (unfinished).
void py_aquad(py::module& pyopensn);
void wrap_quadrature(py::module& aquad);

/// Wrap the field function components of OpenSn (unfinished).
void py_ffunc(py::module& pyopensn);
void wrap_field_function(py::module& ffunc);

/// Wrap the logical volume components of OpenSn (unfinished).
void py_logvol(py::module& pyopensn);
void wrap_logical_volume(py::module& logvol);

/// Wrap the material components of OpenSn.
void py_mat(py::module& pyopensn);
void wrap_material(py::module& mat);

/// Wrap the mesh components of OpenSn (unfinished).
void py_mesh(py::module& pyopensn);
void wrap_mesh(py::module& mesh);
void wrap_mesh_generator(py::module& mesh);
void wrap_graph_partitioner(py::module& mesh);

/// Wrap the solver components of OpenSn (unfinished).
void py_solver(py::module& pyopensn);
void wrap_solver(py::module& slv);

/// Wrap the source components of OpenSn (unfinished).
void py_source(py::module& pyopensn);
void wrap_source(py::module& src);

/// Wrap the cross section components of OpenSn.
void py_xs(py::module& pyopensn);
void wrap_multigroup_xs(py::module& xs);
void wrap_create_load(py::module& xs);

} // namespace opensn
