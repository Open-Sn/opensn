// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/pybind11.h>
#include <vector>

#include "framework/parameters/parameter_block.h"

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

/// Convert a Python sequence to opensn::Vector3
inline Vector3 to_vect3(py::sequence & seq)
{
  switch (seq.size()) {
    case 0 :
      return Vector3();
    case 1 :
      return Vector3(seq[0].cast<double>());
    case 2 :
      return Vector3(seq[0].cast<double>(), seq[1].cast<double>());
    case 3 :
      return Vector3(seq[0].cast<double>(), seq[1].cast<double>(), seq[2].cast<double>());
  }
  throw std::range_error("Expected sequence of size less than 4, got " + std::to_string(seq.size()) + ".");
}

/// Wrap the angular quadrature components of OpenSn (unfinished).
void py_aquad(py::module& pyopensn);
void WrapQuadrature(py::module& aquad);
void WrapProductQuadrature(py::module& aquad);
void WrapCurvilinearQuadrature(py::module& aquad);
void WrapSLDFESQuadrature(py::module& aquad);

/// Wrap the field function components of OpenSn (unfinished).
void py_ffunc(py::module& pyopensn);
void wrap_field_function(py::module& ffunc);

/// Wrap the logical volume components of OpenSn (unfinished).
void py_logvol(py::module& pyopensn);
void wrap_logical_volume(py::module& logvol);

/// Wrap the material components of OpenSn.
// void py_mat(py::module& pyopensn);
// void WrapMaterial(py::module& mat);

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
void WrapMultiGroupXS(py::module& xs);
void WrapCreateLoadXS(py::module& xs);

} // namespace opensn
