/*
 * Created on Sat, February 22
 *
 * Copyright (c) 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
 */
#ifndef PYOPENSN_PYAPI_HPP_
#define PYOPENSN_PYAPI_HPP_

#include <cstddef>  // std::size_t
#include <vector>   // std::vector

#include <mpi.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "framework/parameters/parameter_block.h"   // opensn::ParameterBlock
#include "framework/parameters/input_parameters.h"  // opensn::InputParameters

namespace opensn {

// Environment for running Python
// ------------------------------

/** @brief Environment initializer and finalizer for OpenSn.*/
class PyEnv {
  public:
    /** @brief Constructor.*/
    PyEnv(void);
    /** @brief Destructor.*/
    ~PyEnv(void);
    /** @brief Default environement.*/
    static PyEnv * p_default_env;
};

// Helpers (transform Python types into C++/OpenSn native type)
// ------------------------------------------------------------

/** @brief Convert a C++ vector to a Python memoryview.*/
template <typename T>
py::memoryview convert_vector(const std::vector<T> & vec) {
    return py::memoryview::from_buffer(const_cast<T *>(vec.data()), {vec.size()}, {sizeof(T)}, true);
}

/** @brief Translate a Python dictionary into a ParameterBlock.*/
ParameterBlock kwargs_to_param_block(const py::kwargs & params);

// Wrap libraries
// --------------

/** @brief Wrap the angular quadrature components of OpenSn (unfinished).*/
void py_aquad(py::module & pyopensn);

/** @brief Wrap the field function components of OpenSn (unfinished).*/
void py_ffunc(py::module & pyopensn);

/** @brief Wrap the linear Boltzmann solver components of OpenSn (unfinished).*/
void py_lbs(py::module & pyopensn);

/** @brief Wrap the logical volume components of OpenSn (unfinished).*/
void py_logvol(py::module & pyopensn);

/** @brief Wrap the material components of OpenSn.*/
void py_mat(py::module & pyopensn);

/** @brief Wrap the mesh components of OpenSn (unfinished).*/
void py_mesh(py::module & pyopensn);

/** @brief Wrap the solver components of OpenSn (unfinished).*/
void py_solver(py::module & pyopensn);

/** @brief Wrap the source components of OpenSn (unfinished).*/
void py_source(py::module & pyopensn);

/** @brief Wrap the cross section components of OpenSn.*/
void py_xs(py::module & pyopensn);

}  // namespace opensn

#endif  // PYOPENSN_PYAPI_HPP_
