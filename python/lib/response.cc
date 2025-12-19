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
    lbs_problem: pyopensn.solver.Solver
        A handle to an existing LBS problem.
    options: Dict
        The specification of adjoint buffers and forward sources to use. Supported keys:
          - buffers: List[Dict], default=[]
              Each dictionary supports:
                - name: str (required)
                    Buffer name used for response queries.
                - file_prefixes: Dict (required)
                    Prefix map with optional keys:
                      - flux_moments: str
                          Prefix for flux moments files.
                      - angular_fluxes: str
                          Prefix for angular flux files.
          - clear_sources: bool, default=False
              If true, clears any existing forward sources.
          - sources: Dict, default={}
              Dictionary describing forward sources. Supported keys:
                - material: List[Dict], default=[]
                    Each dictionary supports:
                      - block_id: int (required)
                          Block id for the material source.
                      - strength: List[float] (required)
                          Group-wise source strength (length must match number of groups).
                - point: List[pyopensn.source.PointSource], default=[]
                    Point source handles.
                - volumetric: List[pyopensn.source.VolumetricSource], default=[]
                    Volumetric source handles.
                - boundary: List[Dict], default=[]
                    Boundary source dictionaries. Each dictionary supports:
                      - name: str (required)
                          Boundary name that identifies the specific boundary.
                      - type: {'vacuum', 'isotropic', 'reflecting'} (required)
                          Boundary type specification. Only ``isotropic`` is used for responses.
                      - group_strength: List[float], optional
                          Required when ``type='isotropic'``. Isotropic strength per group.
                      - function_name: str, optional, default=''
                          Name of a registered function to evaluate on the boundary.
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
        Adjoint buffer name
    )",
    py::arg("buffer_name")
  );
  res_eval.def(
    "SetOptions",
    [](ResponseEvaluator& self, py::kwargs& params)
    {
      InputParameters input = ResponseEvaluator::GetOptionsBlock();
      input.AssignParameters(kwargs_to_param_block(params));
      self.SetOptions(input);
    },
    R"(
    Set options for the response evaluator for adding adjoint buffers and defining forward sources.

    Parameters
    ----------
    buffers: List[Dict], default=[]
        A list of dictionaries containing adjoint buffer specifications. Each dictionary supports:
          - name: str (required)
              Buffer name used for response queries.
          - file_prefixes: Dict (required)
              Prefix map with optional keys:
                - flux_moments: str
                    Prefix for flux moments files.
                - angular_fluxes: str
                    Prefix for angular flux files.
    clear_sources: bool, default=False
        A flag to clear existing sources.
    sources: Dict, default={}
        Dictionary describing forward sources. Supported keys:
          - material: List[Dict], default=[]
              Each dictionary supports:
                - block_id: int (required)
                    Block id for the material source.
                - strength: List[float] (required)
                    Group-wise source strength (length must match number of groups).
          - point: List[pyopensn.source.PointSource], default=[]
              Point source handles.
          - volumetric: List[pyopensn.source.VolumetricSource], default=[]
              Volumetric source handles.
          - boundary: List[Dict], default=[]
              Boundary source dictionaries. Each dictionary supports:
                - name: str (required)
                    Boundary name that identifies the specific boundary.
                - type: {'vacuum', 'isotropic', 'reflecting'} (required)
                    Boundary type specification. Only ``isotropic`` is used for responses.
                - group_strength: List[float], optional
                    Required when ``type='isotropic'``. Isotropic strength per group.
                - function_name: str, optional, default=''
                    Name of a registered function to evaluate on the boundary.
    )"
    );
  res_eval.def(
    "ClearForwardSources",
    &ResponseEvaluator::ClearForwardSources,
    "Clear the existing forward sources from the response evaluator."
  );
  res_eval.def(
    "SetBufferOptions",
    [](ResponseEvaluator& self, py::kwargs& params)
    {
      InputParameters input = ResponseEvaluator::GetBufferOptionsBlock();
      input.AssignParameters(kwargs_to_param_block(params));
      self.SetBufferOptions(input);
    },
    R"(
    Add adjoint buffers from the options block to the response evaluator.

    Parameters
    ----------
    name: str (required)
        Buffer name used for response queries.
    file_prefixes: Dict (required)
        Prefix map with optional keys:
          - flux_moments: str
              Prefix for flux moments files.
          - angular_fluxes: str
              Prefix for angular flux files.
    )"
  );
  res_eval.def(
    "SetSourceOptions",
    [](ResponseEvaluator& self, py::kwargs& params)
    {
      InputParameters input = ResponseEvaluator::GetSourceOptionsBlock();
      input.AssignParameters(kwargs_to_param_block(params));
      self.SetSourceOptions(input);
    },
    R"(
    Add forward sources from the options block to the response evaluator.

    Parameters
    ----------
    material: List[Dict], default=[]
        Each dictionary supports:
          - block_id: int (required)
              Block id for the material source.
          - strength: List[float] (required)
              Group-wise source strength (length must match number of groups).
    point: List[pyopensn.source.PointSource], default=[]
        Point source handles.
    volumetric: List[pyopensn.source.VolumetricSource], default=[]
        Volumetric source handles.
    boundary: List[Dict], default=[]
        Boundary source dictionaries. Each dictionary supports:
          - name: str (required)
              Boundary name that identifies the specific boundary.
          - type: {'vacuum', 'isotropic', 'reflecting'} (required)
              Boundary type specification. Only ``isotropic`` is used for responses.
          - group_strength: List[float], optional
              Required when ``type='isotropic'``. Isotropic strength per group.
          - function_name: str, optional, default=''
              Name of a registered function to evaluate on the boundary.
    )"
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
