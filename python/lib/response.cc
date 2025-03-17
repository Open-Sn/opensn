// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/response_evaluator.h"
#include <memory>

namespace opensn
{

// Wrap response evaluator
void
WrapResEval(py::module& response)
{
  // clang-format off
  // response evaluator
  auto res_eval = py::class_<ResponseEvaluator, std::shared_ptr<ResponseEvaluator>>(
    response,
    "ResponseEvaluator",
    R"(
    Response evaluator by folding sources against adjoint solutions.

    Wrapper of :cpp:class:`opensn::ResponseEvaluator`.

    The workflow for this utility is constructed to minimize the file reading necessary for
    evaluations. To begin, one should add all adjoint solutions that are desired for response
    computations into the buffer. Then, one should define the different forward source
    configurations of interest. With this, the user can now iterate over the source
    configurations and convolve them against all available adjoint solutions in the buffer.
    )"
  );
  res_eval.def(
    py::init(
      [](py::kwargs& params)
      {
        return ResponseEvaluator::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a response evaluator object.

    Parameters
    ----------
    ???
    )"
  );
  res_eval.def(
    "EvaluateResponse",
    &ResponseEvaluator::EvaluateResponse,
    R"(
    Evaluate a response using the specified adjoint buffer with the currently defined sources in
    the solver.

    Parameters
    ----------
    buffer_name: str
        ???
    )",
    py::arg("buffer_name")
  );
  res_eval.def(
    "ClearForwardSources",
    &ResponseEvaluator::ClearForwardSources,
    "Clear the existing forward sources from the response evaluator."
  );
  // clang-format on
}

// Wrap the response components of OpenSn
void
py_response(py::module& pyopensn)
{
  py::module response = pyopensn.def_submodule("response", "Response module.");
  WrapResEval(response);
}

} // namespace opensn
