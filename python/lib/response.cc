// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/response_evaluator.h"
#include <pybind11/stl.h>
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
    computations as named data sets. Then, one should define the different forward source
    configurations of interest. With this, the user can now iterate over the source
    configurations and convolve them against all available adjoint data sets.
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
    problem: pyopensn.solver.DiscreteOrdinatesProblem
        A handle to an existing LBS problem.
    options: Dict
        The specification of data sets and forward sources to use. Supported keys:
          - data: List[Dict], default=[]
              Data set specifications, as for :py:meth:`SetOptions`.
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
    Evaluate a response by folding the currently defined forward sources against the adjoint
    flux moments and angular fluxes of a data set added with ``SetOptions(data=...)`` or
    :py:meth:`SetDataOptions`. This call is collective.

    Parameters
    ----------
    data_name: str
        Name of the data set holding the adjoint flux moments (``flux_moments``) and, for boundary
        sources, adjoint angular fluxes (``angular_fluxes``).

    Raises
    ------
    ValueError
        If no data set has this name.
    RuntimeError
        If the data set lacks the flux moments or angular fluxes the defined sources require.
    )",
    py::arg("data_name")
  );
  res_eval.def(
    "EvaluateSurfaceResponse",
    [](const ResponseEvaluator& self,
       const std::string& fwd_data_name,
       const std::string& adj_data_name,
       const std::vector<std::string>& surfaces)
    {
      py::dict result;
      for (const auto& surface_response :
           self.EvaluateSurfaceResponse(fwd_data_name, adj_data_name, surfaces))
      {
        py::dict parts;
        parts["outgoing"] = surface_response.outgoing;
        parts["incoming"] = surface_response.incoming;
        parts["total"] = surface_response.Total();
        result[py::str(surface_response.name)] = std::move(parts);
      }
      return result;
    },
    R"(
    Evaluate the surface inner products of forward and adjoint surface angular fluxes.

    Both data sets must have been added with ``SetOptions(data=...)`` using the
    ``surface_angular_fluxes`` prefix. Faces are matched by location, so the two data sets may have
    been written by different numbers of ranks and partitions, and neither needs to match the
    partitioning of the current run.

    For a surface :math:`S` with faces :math:`f`, the response carried across it is

    .. math::

       R = \sum_{f \in S} \sum_{d \in D_f} \mu_{f,d}\, w_d \sum_{ij} M_{f,ij} \sum_g
           \psi_{f,i,d,g}\, \psi^\dagger_{f,j,d,g},

    with :math:`\mu_{f,d} = \Omega_d \cdot n_f`, :math:`w_d` the quadrature weight,
    :math:`M_{f,ij}` the mass matrix of face :math:`f`, and :math:`\psi_{f,i,d,g}` and
    :math:`\psi^\dagger_{f,j,d,g}` the forward and adjoint angular fluxes at nodes :math:`i` and
    :math:`j` of face :math:`f`.

    For a boundary, :math:`n_f` is the outward normal and :math:`D_f` holds every direction. For the
    outflow directions (:math:`\mu_{f,d} > 0`), the forward flux is that of the face's cell and the
    adjoint flux is the adjoint boundary-condition value (0 for a vacuum boundary), which is upwind
    for the adjoint because it travels along :math:`-\Omega_d`. For the inflow directions, the
    forward flux is the forward boundary-condition value (0 for a vacuum boundary) and the adjoint
    flux that of the face's cell. This is the boundary term of the discrete adjoint identity, so a
    boundary response is exact for the discretization, like an interior surface response: its
    ``outgoing`` part is the response to the particles leaving through the boundary (for a
    detector outside the domain, represented by an adjoint boundary source, the detector
    response), and its ``incoming`` part minus the response to a forward boundary source on it.

    The surface files record whether each tag is a boundary or one side of an interior surface,
    so a boundary is evaluated as a boundary whatever its name, including names ending in ``_u``
    or ``_d``.

    An interior surface is stored under tags ``<name>_u`` and ``<name>_d``, one per side. Its faces
    are those of the ``_u`` side, each matched to the ``_d`` face at the same location, and
    :math:`n_f` is the outward normal of the side the surface is requested by: the ``_u`` normal
    (the positive axis) for ``<name>`` or ``<name>_u``, and the ``_d`` normal (the negative axis)
    for ``<name>_d``, which gives the same response with the opposite sign. :math:`D_f` holds
    every direction. The forward flux is taken from the side of the face upwind of
    :math:`\Omega_d`, the side it leaves. The adjoint flux is taken from the downwind side, which
    is upwind for the adjoint because it travels along :math:`-\Omega_d`.

    When an interior surface separates the forward source from the detector (the region the adjoint
    source scores), its response equals the detector response if :math:`n_f` points from the source
    side to the detector side, and the negative of it otherwise. Request the surface by the tag of
    the side the source is on to get the positive value: ``<name>_u`` when the source is on the
    negative side of the surface, ``<name>_d`` when it is on the positive side. A surface that does
    not separate them gives zero.

    Each response is split by the sign of :math:`\mu_{f,d}`: ``outgoing`` sums the directions
    with :math:`\mu_{f,d} > 0`, which cross along :math:`n_f`, and ``incoming`` those with
    :math:`\mu_{f,d} < 0`, which cross back; ``total`` is their sum, :math:`R`. For an interior
    surface ``x20`` requested along its ``_u`` normal (+x), ``outgoing`` is the part carried in the
    +x direction and ``incoming`` minus the part carried in the -x direction. Requested as
    ``x20_d``, the two parts swap and change sign.

    Forward and adjoint directions are paired by equal direction vectors, and the nodes of matched
    faces by location. The responses are also logged at verbosity 1. This call is collective.
    The adjoint data must come from an adjoint problem, which OpenSn supports only in Cartesian
    geometry, so surface responses are not available for curvilinear (RZ) problems.

    Parameters
    ----------
    fwd_data_name: str
        Name of the forward data set. It must have been written from a forward problem.
    adj_data_name: str
        Name of the adjoint data set. It must have been written from an adjoint problem.
    surfaces: List[str], default=[]
        Surfaces to evaluate: boundary names, and interior surface names, bare (``<name>``, along
        the ``_u`` normal) or with the tag of the side whose normal to use (``<name>_u`` or
        ``<name>_d``). Each surface may be requested once, under one of its names. If empty, every
        surface in the forward data is evaluated, in name order, interior surfaces along their
        ``_u`` normals.

    Returns
    -------
    Dict[str, Dict[str, float]]
        For each surface, under the name it was requested by and in the requested order, a
        dictionary with keys ``outgoing``, ``incoming``, and ``total`` (their sum, the response
        carried across the surface). Sum the ``total`` values to combine surfaces.

    Raises
    ------
    ValueError
        If a data set or surface is unknown, a surface is requested more than once, a required tag
        is missing from the adjoint data or recorded there as a different kind of surface, or the
        data sets have incompatible quadratures, groupsets, adjoint flags, or face keys.
    RuntimeError
        If the forward and adjoint data do not cover the same faces.
    )",
    py::arg("fwd_data_name"),
    py::arg("adj_data_name"),
    py::arg("surfaces") = std::vector<std::string>{}
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
    Set options for the response evaluator for adding data sets and defining forward sources.

    Parameters
    ----------
    data: List[Dict], default=[]
        A list of dictionaries containing data set specifications, each as accepted by
        :py:meth:`SetDataOptions`.
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
    "SetDataOptions",
    [](ResponseEvaluator& self, py::kwargs& params)
    {
      InputParameters input = ResponseEvaluator::GetDataOptionsBlock();
      input.AssignParameters(kwargs_to_param_block(params));
      self.SetDataOptions(input);
    },
    R"(
    Add a named data set to the response evaluator. This call is collective.

    Parameters
    ----------
    name: str (required)
        Data set name used by :py:meth:`EvaluateResponse` and :py:meth:`EvaluateSurfaceResponse`.
        It must not already be in use.
    file_prefixes: Dict (required)
        Prefix map with optional keys:
          - flux_moments: str
              Prefix for the adjoint flux-moment files that :py:meth:`EvaluateResponse` folds
              against the forward sources. Written with the current rank count and partitioning.
          - angular_fluxes: str
              Prefix for the adjoint angular-flux files needed for boundary sources. Written with
              the current rank count and partitioning.
          - surface_angular_fluxes: str
              Prefix for the files written by ``WriteSurfaceAngularFluxes``, used by
              :py:meth:`EvaluateSurfaceResponse`. All files of the prefix are read, whatever
              number of ranks wrote them, and their faces are distributed over the current ranks.
    surfaces: List[str], default=[]
        Surfaces to read. Required with ``surface_angular_fluxes``, and rejected without it. Each
        entry is a boundary name or an interior surface name, which reads both stored sides
        ``<name>_u`` and ``<name>_d``. A single side's tag is rejected, because a surface response
        needs both sides.
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
