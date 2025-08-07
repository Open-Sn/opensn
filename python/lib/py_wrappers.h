// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/parameter_block.h"
#include "framework/data_types/vector.h"
#include <pybind11/pybind11.h>
#include <unordered_set>
#include <unordered_map>
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

// Cast kwargs to ParameterBlock

/// Translate a Python dictionary into a ParameterBlock.
ParameterBlock kwargs_to_param_block(const py::kwargs& params);

// Retrieve arguments from the Python interface

/// Pop an object from kwargs with default value
inline py::object
pop_cast(py::kwargs& kwargs, const std::string& key, const py::object& default_value)
{
  if (kwargs.contains(key.c_str()))
  {
    return kwargs.attr("pop")(key);
  }
  return default_value;
}

/// Pop an object from kwargs and raise error if the key is not found
inline py::object
pop_cast(py::kwargs& kwargs, const std::string& key)
{
  if (!kwargs.contains(key.c_str()))
  {
    throw std::runtime_error("Key \"" + key + "\" must be provided.\n");
  }
  return kwargs.attr("pop")(key);
}

/// Extract tuple of arguments
template <typename... Args>
std::tuple<Args...>
extract_args_tuple(py::kwargs& kwargs,
                   const std::vector<std::string>& required_keys,
                   const std::vector<std::pair<std::string, py::object>>& optional_keys =
                     std::vector<std::pair<std::string, py::object>>())
{
  // check size
  if (required_keys.size() + optional_keys.size() != sizeof...(Args))
  {
    throw std::runtime_error(
      "Mismatch number of arguments. "
      "This is a dev bug, please contact OpenSn developpers for this issue!");
  }
  // initialize retriever
  std::size_t index = 0;
  auto retriever = [&](auto& arg)
  {
    using ArgType = std::decay_t<decltype(arg)>;
    if (index < required_keys.size())
    {
      const std::string& key = required_keys[index];
      arg = pop_cast(kwargs, key).cast<ArgType>();
    }
    else
    {
      const auto& [key, default_val] = optional_keys[index - required_keys.size()];
      arg = pop_cast(kwargs, key, default_val).cast<ArgType>();
    }
    ++index;
  };
  // retrieve keys
  std::tuple<Args...> args;
  std::apply([&](auto&... args) { (retriever(args), ...); }, args);
  // check for orphan keys
  if (!kwargs.empty())
  {
    std::ostringstream err_log;
    err_log << "Unknown arguments(s):";
    for (const auto& item : kwargs)
    {
      std::string key = py::str(item.first);
      err_log << " \"" << key << "\"";
    }
    err_log << ".";
    throw std::runtime_error(err_log.str());
  }
  return args;
}

/**
 *  @brief Construct an object from kwargs
 *  @details This function template allows construction of a C++ object from Python keyword
 *  arguments (`py::kwargs`), with enforcement of required and optional arguments.
 *  @tparam T The target class type to construct.
 *  @tparam Args The constructor argument types of T.
 *  @param kwargs Python keyword arguments provided from a Pybind11 binding.
 *  @param required_keys List of required argument names (must appear in ``kwargs``).
 *  @param optional_keys List of optional arguments with default values (used if not found in
 *  ``kwargs``).
 *  @note This function is meant to replace the functionality of ``InputParameters``.
 *  @example
 *  // C++ class
 *  class Foo {
 *    public:
 *      Foo(int a, double b, const std::string& name, bool verbose = false);
 *  };
 *
 *  // Pybind11 constructor wrapper
 *  foo.def(
 *    py::init(
 *      [](py::kwargs& params) {
 *        const std::vector<std::string> required_keys = {"a", "b", "name"};
 *        const std::vector<std::pair<std::string, py::object>> optional_keys = {
 *          {"verbose", py::bool_(false)}
 *        };
 *        return construct_from_kwargs<Foo, int, double, std::string, bool>(params, required_keys,
 *                                                                          optional_keys);
 *      }
 *    ),
 *    R"(
 *    Construct ...
 *
 *    Parameters
 *    ----------
 *    a: int
 *        ...
 *    b: float
 *        ...
 *    name: str
 *        ...
 *    verbose: bool, default=False
 *        ...
 *    )"
 *  );
 */
template <class T, typename... Args>
std::shared_ptr<T>
construct_from_kwargs(py::kwargs& kwargs,
                      const std::vector<std::string>& required_keys,
                      const std::vector<std::pair<std::string, py::object>>& optional_keys =
                        std::vector<std::pair<std::string, py::object>>())
{
  std::tuple<Args...> args = extract_args_tuple<Args...>(kwargs, required_keys, optional_keys);
  return std::apply(
    [](auto&&... unpacked_args)
    { return std::make_shared<T>(std::forward<decltype(unpacked_args)>(unpacked_args)...); },
    args);
}

// Module wrappers

/// Wrap the context components of OpenSn.
void py_context(py::module& pyopensn);

/// Wrap the angular quadrature components of OpenSn.
void py_aquad(py::module& pyopensn);
void WrapQuadraturePointPhiTheta(py::module& aquad);
void WrapQuadrature(py::module& aquad);
void WrapProductQuadrature(py::module& aquad);
void WrapCurvilinearQuadrature(py::module& aquad);
void WrapSLDFESQuadrature(py::module& aquad);

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

/// Wrap the settings components of OpenSn
void py_settings(py::module& pyopensn);

/// Wrap the solver components of OpenSn (unfinshed).
void py_solver(py::module& pyopensn);
void WrapProblem(py::module& slv);
void WrapSolver(py::module& slv);
void WrapLBS(py::module& slv);
void WrapSteadyState(py::module& slv);
void WrapKEigen(py::module& slv);
void WrapNLKEigen(py::module& slv);
void WrapPIteration(py::module& slv);
void WrapDiscreteOrdinatesKEigenAcceleration(py::module& slv);

/// Wrap the source components of OpenSn.
void py_source(py::module& pyopensn);
void WrapPointSource(py::module& src);
void WrapVolumetricSource(py::module& src);

/// Wrap the cross section components of OpenSn.
void py_xs(py::module& pyopensn);
void WrapMultiGroupXS(py::module& xs);

} // namespace opensn
