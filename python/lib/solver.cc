// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include <pybind11/functional.h>
#include "framework/runtime.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/discrete_ordinates_keigen_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/scdsa_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/smm_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/transient_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/nl_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/solver.h"
#include <pybind11/numpy.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace opensn
{

// Wrap problem
void
WrapProblem(py::module& slv)
{
  // clang-format off
  // problem base
  auto problem = py::class_<Problem, std::shared_ptr<Problem>>(
    slv,
    "Problem",
    R"(
    Base class for all problems.

    Wrapper of :cpp:class:`opensn::Problem`.
    )"
  );
  // clang-format on
}

// Wrap solver
void
WrapSolver(py::module& slv)
{
  // clang-format off
  // solver base
  auto solver = py::class_<Solver, std::shared_ptr<Solver> >(
    slv,
    "Solver",
    R"(
    Base class for all solvers.

    Wrapper of :cpp:class:`opensn::Solver`.
    )"
  );
  solver.def(
    "Initialize",
    &Solver::Initialize,
    "Initialize the solver."
  );
  solver.def(
    "Execute",
    &Solver::Execute,
    "Execute the solver."
  );
  solver.def(
    "Advance",
    &Solver::Advance,
    "Advance time values function."
  );
  // clang-format on
}

// Wrap LBS solver
void
WrapLBS(py::module& slv)
{
  // clang-format off
  // LBS problem
  auto lbs_problem = py::class_<LBSProblem, std::shared_ptr<LBSProblem>, Problem>(
    slv,
    "LBSProblem",
    R"(
    Base class for all linear Boltzmann problems.

    Wrapper of :cpp:class:`opensn::LBSProblem`.
    )"
  );
  lbs_problem.def(
    "GetScalarFluxFieldFunction",
    [](LBSProblem& self, bool only_scalar_flux)
    {
      py::list field_function_list_per_group;
      for (unsigned int group = 0; group < self.GetNumGroups(); group++)
      {
        if (only_scalar_flux)
        {
          field_function_list_per_group.append(self.CreateScalarFluxFieldFunction(group, 0));
        }
        else
        {
          py::list field_function_list_per_moment;
          for (unsigned int moment = 0; moment < self.GetNumMoments(); moment++)
          {
            field_function_list_per_moment.append(self.CreateScalarFluxFieldFunction(group, moment));
          }
          field_function_list_per_group.append(field_function_list_per_moment);
        }
      }
      return field_function_list_per_group;
    },
    R"(
    Return scalar-flux or flux-moment field functions grouped by energy group.

    Parameters
    ----------
    only_scalar_flux : bool, default=True
        If True, returns only the zeroth moment (scalar flux) field function for each group.
        If False, returns all moment field functions for each group.

    Returns
    -------
    Union[List[pyopensn.fieldfunc.FieldFunctionGridBased], List[List[pyopensn.fieldfunc.FieldFunctionGridBased]]]
        If ``only_scalar_flux=True``:
        ``result[g]`` is the scalar-flux field function for group ``g``.

        If ``only_scalar_flux=False``:
        ``result[g][m]`` is the field function for group ``g`` and moment ``m``.

    Notes
    -----
    Field functions are created on demand from the current solver state.

    The returned field functions are snapshots of the solver state at creation time; they are not
    refreshed automatically if the solver state changes.

    They support ``CanUpdate()`` and ``Update()`` while their owning problem is still alive.
    Calling ``Update()`` explicitly refreshes the same field-function object from the solver's
    latest state.

    Calling ``GetScalarFluxFieldFunction(only_scalar_flux=False)`` creates all requested
    moments from the current ``phi`` iterate at the time of the call.

    In the nested form (``only_scalar_flux=False``), the moment index varies fastest
    within each group (inner index = moment, outer index = group).
    )",
    py::arg("only_scalar_flux") = true
  );
  lbs_problem.def(
    "CreateFieldFunction",
    static_cast<std::shared_ptr<FieldFunctionGridBased> (LBSProblem::*)(
      const std::string&, const std::string&, double)>(&LBSProblem::CreateFieldFunction),
    R"(
    Create a named scalar field function derived from a 1D XS or ``power``.

    Parameters
    ----------
    name : str
        Name to assign to the returned field function.
    xs_name : str
        Built-in 1D XS name, custom XS name, or the special value ``power``.
    power_normalization_target : float, default=-1.0
        If positive, scale the derived field function so that the raw power field would
        integrate to this total power.

    Notes
    -----
    The returned field function is created on demand from the current scalar-flux iterate.
    For ordinary XS names this computes ``sum_g xs[g] * phi_g`` at each node.

    The returned field function is a snapshot of the solver state at creation time; it is not
    refreshed automatically if the solver state changes. It supports ``CanUpdate()`` and
    ``Update()`` while its owning problem is still alive. Calling ``Update()`` explicitly
    recomputes the same field function from the solver's latest state.

    If ``xs_name == "power"``, the same power-generation formula used elsewhere by the solver
    is applied on demand.

    If ``power_normalization_target > 0``, the returned field function is scaled using the power
    implied by the current scalar flux. This scaling affects only the returned field function;
    it does not mutate the solver's internal ``phi`` vectors.

    The returned field function is a fresh object created for this call. It is not
    automatically updated by later solves or timesteps.
    )",
    py::arg("name"),
    py::arg("xs_name"),
    py::arg("power_normalization_target") = -1.0
  );
  lbs_problem.def(
    "GetTime",
    &LBSProblem::GetTime,
    R"(
    Get the current simulation time in seconds.
    )"
  );
  lbs_problem.def(
    "GetTimeStep",
    &LBSProblem::GetTimeStep,
    R"(
    Get the current timestep size.
    )"
  );
  lbs_problem.def(
    "ComputeFissionRate",
    [](LBSProblem& self, const std::string& scalar_flux_iterate)
    {
      const std::vector<double>* phi_ptr = nullptr;
      if (scalar_flux_iterate == "old")
      {
        phi_ptr = &self.GetPhiOldLocal();
      }
      else if (scalar_flux_iterate == "new")
      {
        phi_ptr = &self.GetPhiNewLocal();
      }
      else
      {
        throw std::invalid_argument("Unknown scalar_flux_iterate value: \"" + scalar_flux_iterate + "\".");
      }
      return ComputeFissionRate(self, *phi_ptr);
    },
    R"(
    Computes the total fission rate.

    Parameters
    ----------
    scalar_flux_iterate : {'old', 'new'}
        Specifies which scalar flux vector to use in the calculation.
            - 'old': Use the previous scalar flux iterate.
            - 'new': Use the current scalar flux iterate.

    Returns
    -------
    float
        The total fission rate.

    Raises
    ------
    ValueError
        If `scalar_flux_iterate` is not 'old' or 'new'.
    )",
    py::arg("scalar_flux_iterate")
  );
  lbs_problem.def(
    "ComputeFissionProduction",
    [](LBSProblem& self, const std::string& scalar_flux_iterate)
    {
      const std::vector<double>* phi_ptr = nullptr;
      if (scalar_flux_iterate == "old")
      {
        phi_ptr = &self.GetPhiOldLocal();
      }
      else if (scalar_flux_iterate == "new")
      {
        phi_ptr = &self.GetPhiNewLocal();
      }
      else
      {
        throw std::invalid_argument("Unknown scalar_flux_iterate value: \"" + scalar_flux_iterate + "\".");
      }
      return ComputeFissionProduction(self, *phi_ptr);
    },
    R"(
    Computes the total fission production (nu*fission).

    Parameters
    ----------
    scalar_flux_iterate : {'old', 'new'}
        Specifies which scalar flux vector to use in the calculation.
            - 'old': Use the previous scalar flux iterate.
            - 'new': Use the current scalar flux iterate.

    Returns
    -------
    float
        The total fission production.

    Raises
    ------
    ValueError
        If `scalar_flux_iterate` is not 'old' or 'new'.
    )",
    py::arg("scalar_flux_iterate")
  );
  lbs_problem.def(
    "GetPhiOldLocal",
    [](LBSProblem& self)
    {
      return convert_vector(self.GetPhiOldLocal());
    },
    R"(
    Get the previous scalar flux iterate (local vector).

    Returns
    -------
    memoryview
        Memory view of the local old scalar flux vector.
    )"
  );
  lbs_problem.def(
    "GetPhiNewLocal",
    [](LBSProblem& self)
    {
      return convert_vector(self.GetPhiNewLocal());
    },
    R"(
    Get the current scalar flux iterate (local vector).

    Returns
    -------
    memoryview
        Memory view of the local new scalar flux vector.
    )"
  );
  lbs_problem.def(
    "WriteFluxMoments",
    [](LBSProblem& self, const std::string& file_base)
    {
      LBSSolverIO::WriteFluxMoments(self, file_base);
    },
    R"(
    Write flux moments to file.

    Parameters
    ----------
    file_base: str
        File basename.
    )",
    py::arg("file_base")
  );
  lbs_problem.def(
    "CreateAndWriteSourceMoments",
    [](LBSProblem& self, const std::string& file_base)
    {
      std::vector<double> source_moments = self.MakeSourceMomentsFromPhi();
      LBSSolverIO::WriteFluxMoments(self, file_base, source_moments);
    },
    R"(
    Write source moments from latest flux iterate to file.

    Parameters
    ----------
    file_base: str
        File basename.
    )",
    py::arg("file_base")
  );
  lbs_problem.def(
    "ReadFluxMomentsAndMakeSourceMoments",
    [](LBSProblem& self, const std::string& file_base, bool single_file_flag)
    {
      auto ext_src_moments = self.GetExtSrcMomentsLocal();
      LBSSolverIO::ReadFluxMoments(self, file_base, single_file_flag, ext_src_moments);
      self.SetExtSrcMomentsFrom(ext_src_moments);
      log.Log() << "Making source moments from flux file.";
      const auto temp_phi = self.GetPhiOldLocal();
      self.SetPhiOldFrom(self.GetExtSrcMomentsLocal());
      self.SetExtSrcMomentsFrom(self.MakeSourceMomentsFromPhi());
      self.SetPhiOldFrom(temp_phi);
    },
    R"(
    Read flux moments and compute corresponding source moments.

    Parameters
    ----------
    file_base: str
        File basename.
    single_file_flag: bool
        True if all flux moments are in a single file.
    )",
    py::arg("file_base"),
    py::arg("single_file_flag")
  );
  lbs_problem.def(
    "ReadSourceMoments",
    [](LBSProblem& self, const std::string& file_base, bool single_file_flag)
    {
      auto ext_src_moments = self.GetExtSrcMomentsLocal();
      LBSSolverIO::ReadFluxMoments(self, file_base, single_file_flag, ext_src_moments);
      self.SetExtSrcMomentsFrom(ext_src_moments);
    },
    R"(
    Read source moments from file.

    Parameters
    ----------
    file_base: str
        File basename.
    single_file_flag: bool
        True if all source moments are in a single file.
    )",
    py::arg("file_base"),
    py::arg("single_file_flag")
  );
  lbs_problem.def(
    "ReadFluxMoments",
    [](LBSProblem& self, const std::string& file_base, bool single_file_flag)
    {
      LBSSolverIO::ReadFluxMoments(self, file_base, single_file_flag);
    },
    R"(
    Read flux moment data.

    Parameters
    ----------
    file_base: str
        File basename.
    single_file_flag: bool
        True if all flux moments are in a single file.
    )",
    py::arg("file_base"),
    py::arg("single_file_flag")
  );
  lbs_problem.def(
    "WriteAngularFluxes",
    [](DiscreteOrdinatesProblem& self, const std::string& file_base)
    {
      DiscreteOrdinatesProblemIO::WriteAngularFluxes(self, file_base);
    },
    R"(
    Write angular flux data to file.

    Parameters
    ----------
    file_base: str
        File basename.
    )",
    py::arg("file_base")
  );
  lbs_problem.def(
    "ReadAngularFluxes",
    [](DiscreteOrdinatesProblem& self, const std::string& file_base)
    {
      DiscreteOrdinatesProblemIO::ReadAngularFluxes(self, file_base);
    },
    R"(
    Read angular fluxes from file.

    Parameters
    ----------
    file_base: str
        File basename.
    )",
    py::arg("file_base")
  );
  lbs_problem.def(
    "SetPointSources",
    [](LBSProblem& self, py::kwargs& params)
    {
      for (auto [key, value] : params)
      {
        auto c_key = key.cast<std::string>();
        if (c_key == "clear_point_sources")
          self.ClearPointSources();
        else if (c_key == "point_sources")
        {
          auto sources = value.cast<py::list>();
          for (auto source : sources)
          {
            self.AddPointSource(source.cast<std::shared_ptr<PointSource>>());
          }
        }
        else
          throw std::runtime_error("Invalid argument provided to SetPointSources.\n");
      }
    },
    R"(
    Set or clear point sources.

    Parameters
    ----------
    clear_point_sources: bool, default=False
        If true, all current the point sources of the problem are deleted.
    point_sources: List[pyopensn.source.PointSource]
        List of new point sources to be added to the problem.
    )"
  );
  lbs_problem.def(
    "SetVolumetricSources",
    [](LBSProblem& self, py::kwargs& params)
    {
      for (auto [key, value] : params)
      {
        auto c_key = key.cast<std::string>();
        if (c_key == "clear_volumetric_sources")
          self.ClearVolumetricSources();
        else if (c_key == "volumetric_sources")
        {
          auto sources = value.cast<py::list>();
          for (auto source : sources)
          {
            self.AddVolumetricSource(source.cast<std::shared_ptr<VolumetricSource>>());
          }
        }
        else
          throw std::runtime_error("Invalid argument provided to SetVolumetricSources.\n");
      }
    },
    R"(
    Set or clear volumetric sources.

    Parameters
    ----------
    clear_volumetric_sources: bool, default=False
        If true, all current the volumetric sources of the problem are deleted.
    volumetric_sources: List[pyopensn.source.VolumetricSource]
        List of new volumetric sources to be added to the problem.
    )"
  );
  lbs_problem.def(
    "SetXSMap",
    [](LBSProblem& self, py::kwargs& params)
    {
      BlockID2XSMap xs_map;
      for (auto [key, value] : params)
      {
        auto c_key = key.cast<std::string>();
        if (c_key == "xs_map")
        {
          auto xs_entries = value.cast<py::list>();
          for (auto entry : xs_entries)
          {
            InputParameters xs_entry_pars = LBSProblem::GetXSMapEntryBlock();
            xs_entry_pars.AssignParameters(pyobj_to_param_block("", entry.cast<py::dict>()));
            const auto& block_ids =
              xs_entry_pars.GetParam("block_ids").GetVectorValue<unsigned int>();
            auto xs = xs_entry_pars.GetSharedPtrParam<MultiGroupXS>("xs");
            for (const auto& block_id : block_ids)
              xs_map[block_id] = xs;
          }
        }
        else
          throw std::runtime_error("Invalid argument provided to SetXSMap.\n");
      }
      self.SetBlockID2XSMap(xs_map);
    },
    R"(
    Replace the block-id to cross-section map.

    Parameters
    ----------
    xs_map: List[Dict]
        A list of block-id to cross-section mapping dictionaries. Each dictionary supports:
          - block_ids: List[int] (required)
              Mesh block ids to associate with the cross section.
          - xs: pyopensn.xs.MultiGroupXS (required)
              Cross section object.

    Notes
    -----
    The problem is refreshed immediately after replacing the map. Material metadata,
    precursor storage, GPU carriers, and derived solver state owned by the concrete
    problem are rebuilt for the new cross sections.

    If ``options.use_precursors=True``, this flag remains active across XS-map
    changes even when the current map has no precursor-bearing material. In that
    case delayed-neutron source terms are simply inactive until a later map provides
    precursor data again.

    Existing precursor concentrations are remapped by local cell and precursor-family
    index. When the new material for a cell has fewer precursor families, extra old
    families are dropped. When it has more families, newly introduced families are
    initialized to zero. If a cell changes through a material with zero precursors,
    the precursor history for that cell is discarded and any later reintroduced
    precursor families start from zero.

    If any fissionable material in the new map contains delayed-neutron precursor
    data and ``options.use_precursors=True``, all fissionable materials in the map
    must contain precursor data. Non-fissionable materials may have zero precursors.

    Forward/adjoint mode toggles via :meth:`LBSProblem.SetAdjoint` do not change this map.
    The ``MultiGroupXS`` objects themselves are mutable and shared by pointer. If the same
    ``MultiGroupXS`` handle is shared across multiple problems, toggling adjoint mode in one
    problem also changes the transport mode seen by the others.
    )"
  );
  lbs_problem.def(
    "ZeroPhi",
    [](LBSProblem& self)
    {
      self.ZeroPhi();
    },
    R"(
    Zero the scalar-flux vectors (``phi_old`` and ``phi_new``) in-place.
    )"
  );
  lbs_problem.def(
    "SetAdjoint",
    [](LBSProblem& self, bool adjoint)
    {
      self.SetAdjoint(adjoint);
    },
    py::arg("adjoint") = true,
    R"(
    Set forward/adjoint transport mode.

    Parameters
    ----------
    adjoint: bool, default=True
        ``True`` enables adjoint mode and ``False`` enables forward mode.

    Notes
    -----
    This is one of two supported mode-setting paths in Python:
      1. ``options={'adjoint': ...}`` in the problem constructor.
      2. ``SetAdjoint(...)`` (this method).

    If this changes mode, OpenSn performs a full mode-transition reset:
      - Materials are reinitialized in the selected mode.
      - Point and volumetric sources are cleared.
      - Boundary conditions are cleared.
      - Scalar and angular flux vectors are zeroed.

    If this is called with the same mode as the current setting, no reset is performed.

    The block-id to cross-section map is preserved across the transition.
    However, the mode change is applied to the mapped ``MultiGroupXS`` objects themselves.
    If those objects are shared with other problems, they observe the same mode toggle.

    After a mode change, reapply the desired driving terms before solving, typically:
      - :meth:`LBSProblem.SetPointSources`
      - :meth:`LBSProblem.SetVolumetricSources`
      - :meth:`DiscreteOrdinatesProblem.SetBoundaryOptions`

    This routine is intentionally destructive with respect to source/boundary/flux state
    to avoid hidden coupling between forward and adjoint setups.
    )"
  );
  lbs_problem.def(
    "SetForward",
    &LBSProblem::SetForward,
    R"(
    Set forward transport mode.

    Equivalent to ``SetAdjoint(False)``.
    )"
  );
  lbs_problem.def(
    "IsAdjoint",
    &LBSProblem::IsAdjoint,
    R"(
    Return ``True`` if the problem is in adjoint mode, otherwise ``False``.
    )"
  );
  lbs_problem.def(
    "IsTimeDependent",
    &LBSProblem::IsTimeDependent,
    R"(
    Return ``True`` if the problem is in time-dependent mode, otherwise ``False``.
    )"
  );

  // discrete ordinate solver
  auto do_problem = py::class_<DiscreteOrdinatesProblem, std::shared_ptr<DiscreteOrdinatesProblem>,
                               LBSProblem>(
    slv,
    "DiscreteOrdinatesProblem",
    R"(
    Base class for discrete ordinates problems in Cartesian geometry.

    This class implements the algorithms necessary to solve a problem using the discrete ordinates method.
    When paired with a solver base, the result is a solver instance configured for a specific problem type
    (steady-state, transient, adjoint, k-eigenvalue, etc.).

    Wrapper of :cpp:class:`opensn::DiscreteOrdinatesProblem`.
    )"
  );
  do_problem.def(
    py::init(
      [](py::kwargs& params)
      {
        return DiscreteOrdinatesProblem::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a discrete ordinates problem with Cartesian geometry.

    Parameters
    ----------
    mesh : MeshContinuum
        The spatial mesh.
    num_groups : int
        The total number of energy groups.
    groupsets : List[Dict], default=[]
        A list of input parameter blocks, each block provides the iterative properties for a
        groupset. Each dictionary supports:
          - groups_from_to: List[int] (required)
              Two-entry list with the first and last group id for the groupset, e.g. ``[0, 3]``.
          - angular_quadrature: pyopensn.aquad.AngularQuadrature, optional
              Handle to an angular quadrature.
          - angle_aggregation_type: {'polar', 'single', 'azimuthal'}, default='polar'
              Angle aggregation method to use during sweeping.
          - angle_aggregation_num_subsets: int, default=1
              Number of angle subsets used for aggregation.
          - inner_linear_method: {'classic_richardson', 'petsc_richardson',
            'petsc_gmres', 'petsc_bicgstab'}, default='petsc_richardson'
              Iterative method used for inner linear solves.
          - l_abs_tol: float, default=1.0e-6
              Inner linear solver absolute residual tolerance.
          - l_max_its: int, default=200
              Inner linear solver maximum iterations.
          - gmres_restart_interval: int, default=30
              GMRES restart interval, if GMRES is used.
          - allow_cycles: bool, default=True
              Whether cyclic dependencies are allowed in sweeps.
          - apply_wgdsa: bool, default=False
              Enable within-group DSA for this groupset.
          - wgdsa_l_abs_tol: float, default=1.0e-4
              WGDSA linear absolute tolerance.
          - wgdsa_l_max_its: int, default=30
              WGDSA maximum iterations.
          - wgdsa_verbose: bool, default=False
              Verbose WGDSA output.
          - wgdsa_petsc_options: str, default=''
              PETSc options string for the WGDSA solver.
          - apply_tgdsa: bool, default=False
              Enable two-grid DSA for this groupset.
          - tgdsa_l_abs_tol: float, default=1.0e-4
              TGDSA linear absolute tolerance.
          - tgdsa_l_max_its: int, default=30
              TGDSA maximum iterations.
          - tgdsa_verbose: bool, default=False
              Verbose TGDSA output.
          - tgdsa_petsc_options: str, default=''
              PETSc options string for the TGDSA solver.
    xs_map : List[Dict], default=[]
        A list of mappings from block ids to cross-section definitions. Each dictionary supports:
          - block_ids: List[int] (required)
              Mesh block IDs to associate with the cross section.
          - xs: pyopensn.xs.MultiGroupXS (required)
              Cross-section object to assign to the specified blocks.
    boundary_conditions: List[Dict], default=[]
        A list containing tables for each boundary specification. Each dictionary supports:
          - name: str (required)
              Boundary name that identifies the specific boundary.
          - type: {'vacuum', 'isotropic', 'reflecting', 'arbitrary'} (required)
              Boundary type specification.
          - group_strength: List[float], optional
              Required when ``type='isotropic'``. Isotropic strength per group.
          - start_time: float, optional
              Active start time for isotropic boundaries only. Defaults to -infinity.
          - end_time: float, optional
              Active end time for isotropic boundaries only. Defaults to infinity.
          - function: AngularFluxFunction, optional
              Required when ``type='arbitrary'`` unless ``time_function`` is supplied. Callable
              that returns incoming angular flux from group and direction.
          - time_function: AngularFluxTimeFunction, optional
              Required when ``type='arbitrary'`` unless ``function`` is supplied. Callable that
              returns incoming angular flux from group, direction, and time.
        Isotropic boundaries may use ``start_time``/``end_time`` for simple on/off behavior.
        Arbitrary boundaries must specify exactly one of ``function`` or ``time_function``; use
        ``time_function`` for time-dependent arbitrary inflow and handle any active window inside
        that callback.
    point_sources: List[pyopensn.source.PointSource], default=[]
        A list of point sources.
    volumetric_sources: List[pyopensn.source.VolumetricSource], default=[]
        A list of volumetric sources.
    options : Dict, default={}
        A block of optional configuration parameters, including:
          - max_mpi_message_size: int, default=32768
          - restart_writes_enabled: bool, default=False
          - write_delayed_psi_to_restart: bool, default=True
          - read_restart_path: str, default=''
          - write_restart_path: str, default=''
          - write_restart_time_interval: int, default=0
            (must be 0 or >=30)
          - use_precursors: bool, default=True
            Enable delayed-neutron precursor treatment. This is treated as user intent and remains
            active across later ``SetXSMap`` calls, even if the current XS map temporarily contains
            no precursor-bearing material. When XS maps are swapped, existing precursor
            concentrations are remapped by cell and precursor-family index; new families start at
            zero and removed families are discarded. If any fissionable material in the active map
            has precursor data and ``use_precursors=True``, all fissionable materials in that map
            must have precursor data. Non-fissionable materials may have zero precursors.
          - use_source_moments: bool, default=False
          - save_angular_flux: bool, default=False
            Store angular flux state (`psi`) for transient mode, angular-flux
            field functions, and angular-flux I/O.
          - adjoint: bool, default=False
          - verbose_inner_iterations: bool, default=True
            Print inner iteration details, including WGS and AGS iterations.
          - verbose_outer_iterations: bool, default=True
            Print outer solver progress, including PI/NLKE iterations and transient steps.
          - max_ags_iterations: int, default=100
          - ags_tolerance: float, default=1.0e-6
          - ags_convergence_check: {'l2', 'pointwise'}, default='l2'
          - power_default_kappa: float, default=3.20435e-11
          - field_function_prefix_option: {'prefix', 'solver_name'}, default='prefix'
          - field_function_prefix: str, default=''
        These options are applied at problem creation.
    sweep_type : str, default="AAH"
        The sweep type to use. Must be one of `AAH` or `CBC`. Defaults to `AAH`.
    time_dependent : bool, default=False
        If true, the problem starts in time-dependent mode. Otherwise it starts in
        steady-state mode. Requires ``options.save_angular_flux=True``.
    use_gpus : bool, default=False
        A flag specifying whether GPU acceleration is used for the sweep. Currently, only ``AAH`` is
        supported.
    )"
  );
  do_problem.def(
    "SetTimeDependentMode",
    &DiscreteOrdinatesProblem::SetTimeDependentMode,
    R"(
    Set the problem to time-dependent mode.

    Notes
    -----
    Switch problem from steady-state to time-dependent mode. This updates problem
    internals (sweep chunk mode and source-function) while preserving user boundary
    conditions and fixed sources.

    Requires ``options.save_angular_flux=True`` at problem creation.
    )"
  );
  do_problem.def(
    "SetSteadyStateMode",
    &DiscreteOrdinatesProblem::SetSteadyStateMode,
    R"(
    Set the problem to steady-state mode.

    Notes
    -----
    Switch problem from time-dependent to steady-state mode. This updates problem
    internals (sweep chunk mode and source-function) while preserving user boundary
    conditions and fixed sources.
    )"
  );
  do_problem.def(
    "IsTimeDependent",
    &DiscreteOrdinatesProblem::IsTimeDependent,
    R"(
    Return ``True`` if the problem is currently in time-dependent mode.
    )"
  );
  do_problem.def(
    "SetBoundaryOptions",
    [](DiscreteOrdinatesProblem& self, py::kwargs& params)
    {
      for (auto [key, value] : params)
      {
        auto c_key = key.cast<std::string>();
        if (c_key == "clear_boundary_conditions")
          self.ClearBoundaries();
        else if (c_key == "boundary_conditions")
        {
          auto boundaries = value.cast<py::list>();
          for (auto boundary : boundaries)
          {
            InputParameters input = DiscreteOrdinatesProblem::GetBoundaryOptionsBlock();
            input.AssignParameters(pyobj_to_param_block("", boundary.cast<py::dict>()));
            self.SetBoundaryOptions(input);
          }
        }
        else
          throw std::runtime_error("Invalid argument provided to SetBoundaryOptions.\n");
      }
    },
    R"(
    Set or clear boundary conditions.

    Parameters
    ----------
    clear_boundary_conditions: bool, default=False
        If true, all current boundary conditions are deleted.
    boundary_conditions: List[Dict]
        A list of boundary condition dictionaries. Each dictionary supports:
          - name: str (required)
              Boundary name that identifies the specific boundary.
          - type: {'vacuum', 'isotropic', 'reflecting', 'arbitrary'} (required)
              Boundary type specification.
          - group_strength: List[float], optional
              Required when ``type='isotropic'``. Isotropic strength per group.
          - start_time: float, optional
              Active start time for isotropic boundaries only. Defaults to -infinity.
          - end_time: float, optional
              Active end time for isotropic boundaries only. Defaults to infinity.
          - function: AngularFluxFunction, optional
              Required when ``type='arbitrary'`` unless ``time_function`` is supplied. Callable
              that returns incoming angular flux from group and direction.
          - time_function: AngularFluxTimeFunction, optional
              Required when ``type='arbitrary'`` unless ``function`` is supplied. Callable that
              returns incoming angular flux from group, direction, and time.
        Isotropic boundaries may use ``start_time``/``end_time`` for simple on/off behavior.
        Arbitrary boundaries must specify exactly one of ``function`` or ``time_function``; use
        ``time_function`` for time-dependent arbitrary inflow and handle any active window inside
        that callback.

    Notes
    -----
    Mode transitions via :meth:`LBSProblem.SetAdjoint` clear all boundary conditions.
    Reapply boundaries with this method before solving in the new mode.
    )"
  );
  do_problem.def(
    "ZeroPsi",
    [](DiscreteOrdinatesProblem& self)
    {
      self.ZeroPsi();
    },
    R"(
    Zero the angular-flux storage.
    )"
  );
  do_problem.def(
    "GetPsi",
    [](DiscreteOrdinatesProblem& self)
    {
      const auto& psi = self.GetPsiNewLocal();
      py::list psi_list;
      for (const auto& vec : psi)
      {
        auto array = py::array_t<double>(static_cast<py::ssize_t>(vec.size()),
                                         vec.data(),
                                         py::cast(self));
        psi_list.append(array);
      }
      return psi_list;
    },
    R"(
    Return psi as a list of NumPy arrays (float64), using zero-copy views into the
    underlying data.
    )"
  );
  do_problem.def(
    "GetAngularFieldFunctionList",
    [](DiscreteOrdinatesProblem& self, py::list groups, py::list angles)
    {
      std::vector<unsigned int> group_ids;
      std::vector<size_t> angle_ids;
      group_ids.reserve(groups.size());
      angle_ids.reserve(angles.size());

      for (py::handle g : groups)
        group_ids.push_back(g.cast<unsigned int>());
      for (py::handle a : angles)
        angle_ids.push_back(a.cast<size_t>());

      auto ff_list = self.CreateAngularFluxFieldFunctionList(group_ids, angle_ids);
      py::list out;
      for (const auto& ff : ff_list)
        out.append(ff);
      return out;
    },
    R"(
    Create field functions for selected angular flux components.

    Note: You must enable angular flux storage (``save_angular_flux=True``) in
    the problem options. For transient problems this is required. Otherwise
    the returned field functions will remain zero.

    Example
    -------
    .. code-block:: python

       solver.Initialize()
       solver.Execute()
       ang_ff = phys.GetAngularFieldFunctionList(groups=[0], angles=[0])

    For transient/time-dependent runs, each field function is still a snapshot. Either call
    this after each timestep to create a fresh object or keep the returned object and call
    ``Update()`` after each timestep before exporting or interpolating it. Two common patterns
    are:

    1) Use ``TransientSolver.Execute()`` with a post-advance callback:

    .. code-block:: python

       solver = TransientSolver(problem=phys)
       solver.Initialize()
       ang_ff = phys.GetAngularFieldFunctionList(groups=[0], angles=[0])

       def post_advance():
           for ff in ang_ff:
               ff.Update()
           FieldFunctionGridBased.ExportMultipleToPVTU(ang_ff, "angular_flux_t")

       solver.SetPostAdvanceCallback(post_advance)
       solver.Execute()

    2) Use a custom Python loop with ``TransientSolver.Advance()``:

    .. code-block:: python

       solver = TransientSolver(problem=phys)
       solver.Initialize()
       ang_ff = phys.GetAngularFieldFunctionList(groups=[0], angles=[0])
       for _ in range(num_steps):
           solver.Advance()
           for ff in ang_ff:
               ff.Update()
           FieldFunctionGridBased.ExportMultipleToPVTU(ang_ff, "angular_flux_t")

    Parameters
    ----------
    groups : List[int]
        Global group indices to export.
    angles : List[int]
        Angle indices within the groupset quadrature for each group.

    Returns
    -------
    List[pyopensn.fieldfunc.FieldFunctionGridBased]
        Field functions for the requested ``(group, angle)`` pairs. Each returned field function
        is a snapshot, but supports ``CanUpdate()`` and ``Update()`` while its owning problem is
        alive.
    )",
    py::arg("groups"),
    py::arg("angles")
  );
  do_problem.def(
    "ComputeLeakage",
    [](DiscreteOrdinatesProblem& self, py::list bnd_names)
    {
      auto grid = self.GetGrid();
      // get the supported boundaries
      std::map<std::string, std::uint64_t> allowed_bd_names = grid->GetBoundaryNameMap();
      std::map<std::uint64_t, std::string> allowed_bd_ids = grid->GetBoundaryIDMap();
      const auto coord_sys = grid->GetCoordinateSystem();
      const auto mesh_type = grid->GetType();
      const auto dim = grid->GetDimension();
      // get the boundaries to parse, preserving user order
      std::vector<std::uint64_t> bndry_ids;
      if (bnd_names.size() > 1)
      {
        for (py::handle name : bnd_names)
        {
          auto sname = name.cast<std::string>();
          if (coord_sys == CoordinateSystemType::CYLINDRICAL && dim == 2)
          {
            if (sname == "xmin" || sname == "xmax" || sname == "ymin" || sname == "ymax")
              throw std::runtime_error("ComputeLeakage: Boundary name '" + sname +
                                       "' is invalid for cylindrical orthogonal meshes. "
                                       "Use rmin, rmax, zmin, zmax.");

            if (mesh_type == MeshType::ORTHOGONAL)
            {
              if (sname == "rmin") sname = "xmin";
              else if (sname == "rmax") sname = "xmax";
              else if (sname == "zmin") sname = "ymin";
              else if (sname == "zmax") sname = "ymax";
            }
          }
          bndry_ids.push_back(allowed_bd_names.at(sname));
        }
      }
      else
      {
        bndry_ids = self.GetGrid()->GetUniqueBoundaryIDs();
      }
      // compute the leakage
      std::map<std::uint64_t, std::vector<double>> leakage = ComputeLeakage(self, bndry_ids);
      const bool rz_ortho = (coord_sys == CoordinateSystemType::CYLINDRICAL &&
                             mesh_type == MeshType::ORTHOGONAL && dim == 2);

      auto to_rz_name = [](const std::string& name)
      {
        if (name == "xmin") return std::string("rmin");
        if (name == "xmax") return std::string("rmax");
        if (name == "ymin") return std::string("zmin");
        if (name == "ymax") return std::string("zmax");
        return name;
      };

      // convert result to native Python
      py::dict result;
      for (const auto& bndry_id : bndry_ids)
      {
        const auto it = leakage.find(bndry_id);
        if (it == leakage.end())
          continue;
        // construct numpy array and copy contents
        const auto& grp_wise_leakage = it->second;
        py::array_t<double> np_vector(py::ssize_t(grp_wise_leakage.size()));
        auto buffer = np_vector.request();
        auto *np_vector_data = static_cast<double*>(buffer.ptr);
        std::copy(grp_wise_leakage.begin(), grp_wise_leakage.end(), np_vector_data);
        std::string name = allowed_bd_ids.at(bndry_id);
        if (rz_ortho)
          name = to_rz_name(name);
        result[py::str(name)] = std::move(np_vector);
      }

      return result;
    },
    R"(
    Compute leakage for the problem.

    Parameters
    ----------
    bnd_names : List[str]
        A list of boundary names for which leakage should be computed.

    Returns
    -------
    Dict[str, numpy.ndarray]
        A dictionary mapping boundary names to group-wise leakage vectors.
        Each array contains the outgoing angular flux (per group) integrated over
        the corresponding boundary surface.

    Raises
    ------
    RuntimeError
        If `save_angular_flux` option was not enabled during problem setup.

    ValueError
        If one or more boundary ids are not present on the current mesh.
    )",
    py::arg("bnd_names")
  );

  // discrete ordinates curvilinear problem
  auto do_curvilinear_problem = py::class_<DiscreteOrdinatesCurvilinearProblem,
                                           std::shared_ptr<DiscreteOrdinatesCurvilinearProblem>,
                                           DiscreteOrdinatesProblem>(
    slv,
    "DiscreteOrdinatesCurvilinearProblem",
    R"(
    Base class for discrete ordinates problems in curvilinear geometry.

    Wrapper of :cpp:class:`opensn::DiscreteOrdinatesCurvilinearProblem`.
    )"
  );
  do_curvilinear_problem.def(
    py::init(
      [](py::kwargs& params)
      {
        return DiscreteOrdinatesCurvilinearProblem::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a discrete ordinates problem for curvilinear geometry.

    Warnings
    --------
       DiscreteOrdinatesCurvilinearProblem is **experimental** and should be used with caution!

    Parameters
    ----------
    mesh : MeshContinuum
        The spatial mesh.
    coord_system : int
        Coordinate system to use. Must be set to 2 (cylindrical coordinates).
    num_groups : int
        The total number of energy groups.
    groupsets : list of dict
        A list of input parameter blocks, each block provides the iterative properties for a
        groupset. Each dictionary supports:
          - groups_from_to: List[int] (required)
              Two-entry list with the first and last group id for the groupset, e.g. ``[0, 3]``.
          - angular_quadrature: pyopensn.aquad.AngularQuadrature, optional
              Handle to an angular quadrature.
          - angle_aggregation_type: {'polar', 'single', 'azimuthal'}, default='polar'
              Angle aggregation method to use during sweeping.
          - angle_aggregation_num_subsets: int, default=1
              Number of angle subsets used for aggregation.
          - inner_linear_method: {'classic_richardson', 'petsc_richardson',
            'petsc_gmres', 'petsc_bicgstab'}, default='petsc_richardson'
              Iterative method used for inner linear solves.
          - l_abs_tol: float, default=1.0e-6
              Inner linear solver absolute residual tolerance.
          - l_max_its: int, default=200
              Inner linear solver maximum iterations.
          - gmres_restart_interval: int, default=30
              GMRES restart interval, if GMRES is used.
          - allow_cycles: bool, default=True
              Whether cyclic dependencies are allowed in sweeps.
          - apply_wgdsa: bool, default=False
              Enable within-group DSA for this groupset.
          - wgdsa_l_abs_tol: float, default=1.0e-4
              WGDSA linear absolute tolerance.
          - wgdsa_l_max_its: int, default=30
              WGDSA maximum iterations.
          - wgdsa_verbose: bool, default=False
              Verbose WGDSA output.
          - wgdsa_petsc_options: str, default=''
              PETSc options string for the WGDSA solver.
          - apply_tgdsa: bool, default=False
              Enable two-grid DSA for this groupset.
          - tgdsa_l_abs_tol: float, default=1.0e-4
              TGDSA linear absolute tolerance.
          - tgdsa_l_max_its: int, default=30
              TGDSA maximum iterations.
          - tgdsa_verbose: bool, default=False
              Verbose TGDSA output.
          - tgdsa_petsc_options: str, default=''
    xs_map : list of dict
        A list of mappings from block ids to cross-section definitions. Each dictionary supports:
          - block_ids: List[int] (required)
              Mesh block IDs to associate with the cross section.
          - xs: pyopensn.xs.MultiGroupXS (required)
              Cross-section object to assign to the specified blocks.
    boundary_conditions: List[Dict], default=[]
        A list containing tables for each boundary specification. Each dictionary supports:
          - name: str (required)
              Boundary name that identifies the specific boundary.
          - type: {'vacuum', 'isotropic', 'reflecting', 'arbitrary'} (required)
              Boundary type specification.
          - group_strength: List[float], optional
              Required when ``type='isotropic'``. Isotropic strength per group.
          - start_time: float, optional
              Active start time for isotropic boundaries only. Defaults to -infinity.
          - end_time: float, optional
              Active end time for isotropic boundaries only. Defaults to infinity.
          - function: AngularFluxFunction, optional
              Required when ``type='arbitrary'`` unless ``time_function`` is supplied. Callable
              that returns incoming angular flux from group and direction.
          - time_function: AngularFluxTimeFunction, optional
              Required when ``type='arbitrary'`` unless ``function`` is supplied. Callable that
              returns incoming angular flux from group, direction, and time.
        Isotropic boundaries may use ``start_time``/``end_time`` for simple on/off behavior.
        Arbitrary boundaries must specify exactly one of ``function`` or ``time_function``; use
        ``time_function`` for time-dependent arbitrary inflow and handle any active window inside
        that callback.
    point_sources: List[pyopensn.source.PointSource], default=[]
        A list of point sources.
    volumetric_sources: List[pyopensn.source.VolumetricSource], default=[]
        A list of volumetric sources.
    options : dict, optional
        A block of optional configuration parameters applied at problem creation, including:
          - max_mpi_message_size: int, default=32768
          - restart_writes_enabled: bool, default=False
          - write_delayed_psi_to_restart: bool, default=True
          - read_restart_path: str, default=''
          - write_restart_path: str, default=''
          - write_restart_time_interval: int, default=0
          - use_precursors: bool, default=True
            Enable delayed-neutron precursor treatment. This remains active across later
            ``SetXSMap`` calls, even if the current map temporarily has no precursor-bearing
            material. When XS maps change, precursor concentrations are remapped by cell and
            precursor-family index; newly introduced families start at zero and removed families
            are discarded.
          - use_source_moments: bool, default=False
          - save_angular_flux: bool, default=False
            Store angular flux state (`psi`) for transient mode, angular-flux
            field functions, and angular-flux I/O.
          - verbose_inner_iterations: bool, default=True
            Print inner iteration details, including WGS and AGS iterations.
          - verbose_outer_iterations: bool, default=True
            Print outer solver progress, including PI/NLKE iterations and transient steps.
          - max_ags_iterations: int, default=100
          - ags_tolerance: float, default=1.0e-6
          - ags_convergence_check: {'l2', 'pointwise'}, default='l2'
          - power_default_kappa: float, default=3.20435e-11
          - field_function_prefix_option: {'prefix', 'solver_name'}, default='prefix'
          - field_function_prefix: str, default=''
    sweep_type : str, optional
        The sweep type to use. Must be one of `AAH` or `CBC`. Defaults to `AAH`.
        If ``time_dependent=True``, ``options.save_angular_flux=True`` is required.
    )"
  );
}

// Wrap steady-state solver
void
WrapSteadyState(py::module& slv)
{
  const auto balance_table_to_dict = [](const BalanceTable& table)
  {
    py::dict values;
    values["absorption_rate"] = table.absorption_rate;
    values["production_rate"] = table.production_rate;
    values["inflow_rate"] = table.inflow_rate;
    values["outflow_rate"] = table.outflow_rate;
    values["balance"] = table.balance;
    if (table.initial_inventory.has_value())
      values["initial_inventory"] = table.initial_inventory.value();
    if (table.final_inventory.has_value())
      values["final_inventory"] = table.final_inventory.value();
    if (table.predicted_inventory_change.has_value())
      values["predicted_inventory_change"] = table.predicted_inventory_change.value();
    if (table.actual_inventory_change.has_value())
      values["actual_inventory_change"] = table.actual_inventory_change.value();
    if (table.inventory_residual.has_value())
      values["inventory_residual"] = table.inventory_residual.value();
    return values;
  };

  // clang-format off
  // steady state solver
  auto steady_state_solver = py::class_<SteadyStateSourceSolver, std::shared_ptr<SteadyStateSourceSolver>,
                                        Solver>(
    slv,
    "SteadyStateSourceSolver",
    R"(
    Steady state solver.

    Wrapper of :cpp:class:`opensn::SteadyStateSourceSolver`.
    )"
  );
  steady_state_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return SteadyStateSourceSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a steady state solver.

    Parameters
    ----------
    problem : pyopensn.solver.LBSProblem
        Existing LBSProblem instance.

    Notes
    -----
    If ``problem.options.read_restart_path`` is set, restart data is read
    during :meth:`Initialize`. If ``problem.options.restart_writes_enabled`` is
    true, a restart dump is written after :meth:`Execute` completes.
    problem : pyopensn.solver.DiscreteOrdinatesProblem
        Existing discrete ordinates problem instance.
    )"
  );
  steady_state_solver.def(
    "ComputeBalanceTable",
    [balance_table_to_dict](const SteadyStateSourceSolver& self)
    {
      return balance_table_to_dict(self.ComputeBalanceTable());
    },
    R"(
    Compute and return the global balance table using the solver's normalization.
    This is a collective operation and must be called on all ranks.

    Returns
    -------
    dict
        Dictionary with the following entries:

        - ``absorption_rate``:
          Global absorption rate, approximately ``integral sigma_a * phi dV`` summed over
          groups and the full domain.
        - ``production_rate``:
          Global volumetric production/source rate used by the solver,
          approximately ``integral Q dV`` summed over groups and the full domain.
        - ``inflow_rate``:
          Global incoming boundary contribution integrated over incoming
          angular flux on boundaries.
        - ``outflow_rate``:
          Global outgoing boundary contribution accumulated from face outflow
          tallies.
        - ``balance``:
          Rate balance,
          ``production_rate + inflow_rate - absorption_rate - outflow_rate``.

    Notes
    -----
    This solver applies no extra normalization to the balance table.
    )"
  );
  // clang-format on
}

// Wrap transient solver
void
WrapTransient(py::module& slv)
{
  const auto balance_table_to_dict = [](const BalanceTable& table)
  {
    py::dict values;
    values["absorption_rate"] = table.absorption_rate;
    values["production_rate"] = table.production_rate;
    values["inflow_rate"] = table.inflow_rate;
    values["outflow_rate"] = table.outflow_rate;
    values["balance"] = table.balance;
    if (table.initial_inventory.has_value())
      values["initial_inventory"] = table.initial_inventory.value();
    if (table.final_inventory.has_value())
      values["final_inventory"] = table.final_inventory.value();
    if (table.predicted_inventory_change.has_value())
      values["predicted_inventory_change"] = table.predicted_inventory_change.value();
    if (table.actual_inventory_change.has_value())
      values["actual_inventory_change"] = table.actual_inventory_change.value();
    if (table.inventory_residual.has_value())
      values["inventory_residual"] = table.inventory_residual.value();
    return values;
  };
  // clang-format off
  auto transient_solver =
    py::class_<TransientSolver, std::shared_ptr<TransientSolver>, Solver>(
      slv,
      "TransientSolver",
      R"(
      Transient solver.

      Wrapper of :cpp:class:`opensn::TransientSolver`.
      )"
    );
  transient_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return TransientSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a transient solver.

    Parameters
    ----------
    pyopensn.solver.DiscreteOrdinatesProblem : DiscreteOrdinatesProblem
        Existing discrete ordinates problem instance.
    dt : float, optional, default=2.0e-3
        Time step size used during the simulation.
    stop_time : float, optional, default=0.1
        Simulation end time.
    theta : float, optional, default=0.5
        Time differencing scheme parameter.
    initial_state : str, optional, default="existing"
        Initial state for the transient solve. Allowed values: existing, zero.
        In "zero" mode, scalar flux vectors are reset to zero.
    verbose : bool, optional, default=True
        Enable verbose logging.

    Notes
    -----
    The associated problem must have ``save_angular_flux=True`` enabled. This
    is required for transient problems.

    If ``problem.options.read_restart_path`` is set, restart data is read
    during :meth:`Initialize`. If ``problem.options.restart_writes_enabled`` is
    true, timed restart dumps may be written during :meth:`Execute` and a final
    restart dump is written when execution completes.
    )"
  );
  transient_solver.def(
    "SetTimeStep",
    &TransientSolver::SetTimeStep,
    R"(
    Set the timestep size used by :meth:`Advance`.

    Parameters
    ----------
    dt : float
        New timestep size.
    )");
  transient_solver.def(
    "SetTheta",
    &TransientSolver::SetTheta,
    R"(
    Set the theta parameter used by :meth:`Advance`.

    Parameters
    ----------
    theta : float
        Theta value between 1.0e-16 and 1.
    )");
  transient_solver.def(
    "Advance",
    &TransientSolver::Advance,
    R"(
    Advance the solver by a single timestep.

    Notes
    -----
    You must call :meth:`Initialize` before calling :meth:`Advance` or
    :meth:`Execute`.
    )");
  transient_solver.def(
    "SetPreAdvanceCallback",
    static_cast<void (TransientSolver::*)(std::function<void()>)>(
      &TransientSolver::SetPreAdvanceCallback),
    R"(
    Register a callback that runs before each advance within :meth:`Execute`.

    Parameters
    ----------
    callback : Optional[Callable[[], None]]
        Function invoked before the solver advances a timestep. Pass None to clear.
        If the callback modifies the timestep, the new value is used for the
        upcoming step.
    )");
  transient_solver.def(
    "SetPreAdvanceCallback",
    static_cast<void (TransientSolver::*)(std::nullptr_t)>(
      &TransientSolver::SetPreAdvanceCallback),
    "Clear the PreAdvance callback by passing None.");
  transient_solver.def(
    "SetPostAdvanceCallback",
    static_cast<void (TransientSolver::*)(std::function<void()>)>(
      &TransientSolver::SetPostAdvanceCallback),
    R"(
    Register a callback that runs after each advance within :meth:`Execute`.

    Parameters
    ----------
    callback : Optional[Callable[[], None]]
        Function invoked after the solver advances a timestep. Pass None to clear.
    )");
  transient_solver.def(
    "SetPostAdvanceCallback",
    static_cast<void (TransientSolver::*)(std::nullptr_t)>(
      &TransientSolver::SetPostAdvanceCallback),
    "Clear the PostAdvance callback by passing None.");
  transient_solver.def(
    "ComputeBalanceTable",
    [balance_table_to_dict](const TransientSolver& self)
    {
      return balance_table_to_dict(self.ComputeBalanceTable());
    },
    R"(
    Compute and return the global balance table using the solver's normalization.
    This is a collective operation and must be called on all ranks.

    Returns
    -------
    dict
        Dictionary with the following entries:

        - ``absorption_rate``:
          Global absorption rate, approximately ``integral sigma_a * phi dV`` summed over
          groups and the full domain.
        - ``production_rate``:
          Global volumetric production/source rate used by the solver,
          approximately ``integral Q dV`` summed over groups and the full domain.
        - ``inflow_rate``:
          Global incoming boundary contribution integrated over incoming
          angular flux on boundaries.
        - ``outflow_rate``:
          Global outgoing boundary contribution accumulated from face outflow
          tallies.
        - ``balance``:
          Rate balance,
          ``production_rate + inflow_rate - absorption_rate - outflow_rate``.
        - ``initial_inventory``:
          Total particle inventory at the start of the timestep, computed as
          ``integral (1 / v_g) * phi_old dV`` summed over groups and the full domain.
        - ``final_inventory``:
          Total particle inventory at the end of the timestep, computed as
          ``integral (1 / v_g) * phi_new dV`` summed over groups and the full domain.
        - ``predicted_inventory_change``:
          Inventory change predicted by the current timestep balance, computed as 
          ``dt * balance``.
        - ``actual_inventory_change``:
          Measured change in total particle inventory over the timestep, computed as
          ``final_inventory - initial_inventory``.
        - ``inventory_residual``:
          Mismatch between the measured and predicted timestep inventory
          changes, computed as
          ``actual_inventory_change - predicted_inventory_change``.

    Notes
    -----
    This solver applies no extra normalization to the balance table.

    The transient inventory terms currently use the end-of-step rate balance to
    estimate the timestep inventory change.
    )"
  );
  slv.attr("BackwardEuler") = 1.0;
  slv.attr("CrankNicolson") = 0.5;
  // clang-format on
}

// Wrap non-linear k-eigen solver
void
WrapNLKEigen(py::module& slv)
{
  const auto balance_table_to_dict = [](const BalanceTable& table)
  {
    py::dict values;
    values["absorption_rate"] = table.absorption_rate;
    values["production_rate"] = table.production_rate;
    values["inflow_rate"] = table.inflow_rate;
    values["outflow_rate"] = table.outflow_rate;
    values["balance"] = table.balance;
    if (table.initial_inventory.has_value())
      values["initial_inventory"] = table.initial_inventory.value();
    if (table.final_inventory.has_value())
      values["final_inventory"] = table.final_inventory.value();
    if (table.predicted_inventory_change.has_value())
      values["predicted_inventory_change"] = table.predicted_inventory_change.value();
    if (table.actual_inventory_change.has_value())
      values["actual_inventory_change"] = table.actual_inventory_change.value();
    if (table.inventory_residual.has_value())
      values["inventory_residual"] = table.inventory_residual.value();
    return values;
  };
  // clang-format off
  // non-linear k-eigen solver
  auto non_linear_k_eigen_solver = py::class_<NonLinearKEigenSolver, std::shared_ptr<NonLinearKEigenSolver>,
                                              Solver>(
    slv,
    "NonLinearKEigenSolver",
    R"(
    Non-linear k-eigenvalue solver.

    Wrapper of :cpp:class:`opensn::NonLinearKEigenSolver`.
    )"
  );
  non_linear_k_eigen_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return NonLinearKEigenSolver::Create(kwargs_to_param_block(params));
      }
        ),
    R"(
    Construct a non-linear k-eigenvalue solver.

    Parameters
    ----------
    problem: pyopensn.solver.DiscreteOrdinatesProblem
        Existing discrete ordinates problem instance.
    nl_abs_tol: float, default=1.0e-8
        Non-linear absolute tolerance.
    nl_rel_tol: float, default=1.0e-8
        Non-linear relative tolerance.
    nl_sol_tol: float, default=1.0e-50
        Non-linear solution tolerance.
    nl_max_its: int, default=50
        Non-linear algorithm maximum iterations.
    l_abs_tol: float, default=1.0e-8
        Linear absolute tolerance.
    l_rel_tol: float, default=1.0e-8
        Linear relative tolerance.
    l_div_tol: float, default=1.0e6
        Linear divergence tolerance.
    l_max_its: int, default=50
        Linear algorithm maximum iterations.
    l_gmres_restart_intvl: int, default=30
        GMRES restart interval.
    l_gmres_breakdown_tol: float, default=1.0e6
        GMRES breakdown tolerance.
    reset_phi0: bool, default=True
        If true, reinitializes scalar fluxes to 1.0.
    num_initial_power_iterations: int, default=0
        Number of initial power iterations before the non-linear solve.
    )"
  );
  non_linear_k_eigen_solver.def(
    "GetEigenvalue",
    &NonLinearKEigenSolver::GetEigenvalue,
    R"(
    Return the current k-eigenvalue.
    )"
  );
  non_linear_k_eigen_solver.def(
    "ComputeBalanceTable",
    [balance_table_to_dict](const NonLinearKEigenSolver& self)
    {
      return balance_table_to_dict(self.ComputeBalanceTable());
    },
    R"(
    Compute and return the global balance table using the solver's normalization.
    This is a collective operation and must be called on all ranks.

    Returns
    -------
    dict
        Dictionary with the following entries:

        - ``absorption_rate``:
          Global absorption rate, approximately ``integral sigma_a * phi dV`` summed over
          groups and the full domain.
        - ``production_rate``:
          Global volumetric production/source rate used by the solver,
          approximately ``integral Q dV`` summed over groups and the full domain.
        - ``inflow_rate``:
          Global incoming boundary contribution integrated over incoming
          angular flux on boundaries.
        - ``outflow_rate``:
          Global outgoing boundary contribution accumulated from face outflow
          tallies.
        - ``balance``:
          Rate balance,
          ``production_rate + inflow_rate - absorption_rate - outflow_rate``.

    Notes
    -----
    For k-eigenvalue balance reporting, this solver scales the production term by
    ``1 / k_eff`` before forming both ``production_rate`` and ``balance``.
    )"
  );
  // clang-format on
}

// Wrap power iteration solvers
void
WrapPIteration(py::module& slv)
{
  const auto balance_table_to_dict = [](const BalanceTable& table)
  {
    py::dict values;
    values["absorption_rate"] = table.absorption_rate;
    values["production_rate"] = table.production_rate;
    values["inflow_rate"] = table.inflow_rate;
    values["outflow_rate"] = table.outflow_rate;
    values["balance"] = table.balance;
    if (table.initial_inventory.has_value())
      values["initial_inventory"] = table.initial_inventory.value();
    if (table.final_inventory.has_value())
      values["final_inventory"] = table.final_inventory.value();
    if (table.predicted_inventory_change.has_value())
      values["predicted_inventory_change"] = table.predicted_inventory_change.value();
    if (table.actual_inventory_change.has_value())
      values["actual_inventory_change"] = table.actual_inventory_change.value();
    if (table.inventory_residual.has_value())
      values["inventory_residual"] = table.inventory_residual.value();
    return values;
  };
  // clang-format off
  // power iteration k-eigen solver
  auto pi_k_eigen_solver = py::class_<PowerIterationKEigenSolver, std::shared_ptr<PowerIterationKEigenSolver>,
                                      Solver>(
    slv,
    "PowerIterationKEigenSolver",
    R"(
    Power iteration k-eigenvalue solver.

    Wrapper of :cpp:class:`opensn::PowerIterationKEigenSolver`.
    )"
  );
  pi_k_eigen_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PowerIterationKEigenSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a power iteration k-eigen solver.

    Parameters
    ----------
    problem: pyopensn.solver.DiscreteOrdinatesProblem
        Existing DiscreteOrdinatesProblem instance.
    acceleration: pyopensn.solver.DiscreteOrdinatesKEigenAcceleration
        Optional DiscreteOrdinatesKEigenAcceleration instance for acceleration.
    max_iters: int, default = 1000
        Maximum power iterations allowed.
    k_tol: float, default = 1.0e-10
        Tolerance on the k-eigenvalue.
    reset_solution: bool, default=True
        If true, initialize flux moments to 1.0.
    reset_phi0: bool, default=True
        If true, reinitializes scalar fluxes to 1.0.

    Notes
    -----
    If ``problem.options.read_restart_path`` is set, restart data is read
    during :meth:`Initialize`. If ``problem.options.restart_writes_enabled`` is
    true, timed restart dumps may be written during the outer iteration loop
    and a final restart dump is written when execution completes.
    )"
  );
  pi_k_eigen_solver.def(
    "GetEigenvalue",
    &PowerIterationKEigenSolver::GetEigenvalue,
    R"(
    Return the current k-eigenvalue.
    )"
  );
  pi_k_eigen_solver.def(
    "ComputeBalanceTable",
    [balance_table_to_dict](const PowerIterationKEigenSolver& self)
    {
      return balance_table_to_dict(self.ComputeBalanceTable());
    },
    R"(
    Compute and return the global balance table using the solver's normalization.
    This is a collective operation and must be called on all ranks.

    Returns
    -------
    dict
        Dictionary with the following entries:

        - ``absorption_rate``:
          Global absorption rate, approximately ``integral sigma_a * phi dV`` summed over
          groups and the full domain.
        - ``production_rate``:
          Global volumetric production/source rate used by the solver,
          approximately ``integral Q dV`` summed over groups and the full domain.
        - ``inflow_rate``:
          Global incoming boundary contribution integrated over incoming
          angular flux on boundaries.
        - ``outflow_rate``:
          Global outgoing boundary contribution accumulated from face outflow
          tallies.
        - ``balance``:
          Rate balance,
          ``production_rate + inflow_rate - absorption_rate - outflow_rate``.

    Notes
    -----
    For k-eigenvalue balance reporting, this solver scales the production term by
    ``1 / k_eff`` before forming both ``production_rate`` and ``balance``.
    )"
  );
  // clang-format on
}

// Wrap LBS solver
void
WrapDiscreteOrdinatesKEigenAcceleration(py::module& slv)
{
  // clang-format off
  // discrete ordinates k-eigen acceleration base
  auto acceleration = py::class_<DiscreteOrdinatesKEigenAcceleration,
                                 std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>>(
    slv,
    "DiscreteOrdinatesKEigenAcceleration",
    R"(
    Base class for discrete ordinates k-eigenvalue acceleration methods.

    Wrapper of :cpp:class:`opensn::DiscreteOrdinatesKEigenAcceleration`.
    )"
  );
  // SCDSA acceleration
  auto scdsa_acceleration = py::class_<SCDSAAcceleration,
                                       std::shared_ptr<SCDSAAcceleration>,
                                       DiscreteOrdinatesKEigenAcceleration>(
    slv,
    "SCDSAAcceleration",
    R"(
    Construct an SCDSA accelerator for the power iteration k-eigenvalue solver.

    Wrapper of :cpp:class:`opensn::SCDSAAcceleration`.
    )"
  );
  scdsa_acceleration.def(
    py::init(
      [](py::kwargs& params)
      {
        return SCDSAAcceleration::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    SCDSA acceleration for the power iteration k-eigenvalue solver.

    Parameters
    ----------
    problem: pyopensn.solver.DiscreteOrdinatesProblem
        Existing DiscreteOrdinatesProblem instance.
    l_abs_tol: float, defauilt=1.0e-10
        Absolute residual tolerance.
    max_iters: int, default=100
        Maximum allowable iterations.
    verbose: bool, default=False
        If true, enables verbose output.
    petsc_options: str, default="ssss"
        Additional PETSc options.
    pi_max_its: int, default=50
        Maximum allowable iterations for inner power iterations.
    pi_k_tol: float, default=1.0e-10
        k-eigenvalue tolerance for the inner power iterations.
    sdm: str, default="pwld"
        Spatial discretization method to use for the diffusion solver. Valid choices are:
            - 'pwld' : Piecewise Linear Discontinuous
            - 'pwlc' : Piecewise Linear Continuous
    )"
  );
  // SMM acceleration
  auto smm_acceleration = py::class_<SMMAcceleration,
                                     std::shared_ptr<SMMAcceleration>,
                                     DiscreteOrdinatesKEigenAcceleration>(
    slv,
    "SMMAcceleration",
    R"(
    Construct an SMM accelerator for the power iteration k-eigenvalue solver.

    Wrapper of :cpp:class:`opensn::SMMAcceleration`.
    )"
  );
  smm_acceleration.def(
    py::init(
      [](py::kwargs& params)
      {
        return SMMAcceleration::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    SMM acceleration for the power iteration k-eigenvalue solver.

    Warnings
    --------
       SMM acceleration is **experimental** and should be used with caution!
       SMM accleration only supports problems with isotropic scattering.

    Parameters
    ----------
    problem: pyopensn.solver.DiscreteOrdinatesProblem
        Existing DiscreteOrdinatesProblem instance.
    l_abs_tol: float, defauilt=1.0e-10
        Absolute residual tolerance.
    max_iters: int, default=100
        Maximum allowable iterations.
    verbose: bool, default=False
        If true, enables verbose output.
    petsc_options: str, default="ssss"
        Additional PETSc options.
    pi_max_its: int, default=50
        Maximum allowable iterations for inner power iterations.
    pi_k_tol: float, default=1.0e-10
        k-eigenvalue tolerance for the inner power iterations.
    sdm: str, default="pwld"
        Spatial discretization method to use for the diffusion solver. Valid choices are:
            - 'pwld' : Piecewise Linear Discontinuous
            - 'pwlc' : Piecewise Linear Continuous
    )"
  );
  // clang-format on
}

// Wrap the solver components of OpenSn
void
py_solver(py::module& pyopensn)
{
  py::module slv = pyopensn.def_submodule("solver", "Solver module.");
  WrapProblem(slv);
  WrapSolver(slv);
  WrapLBS(slv);
  WrapSteadyState(slv);
  WrapTransient(slv);
  WrapNLKEigen(slv);
  WrapDiscreteOrdinatesKEigenAcceleration(slv);
  WrapPIteration(slv);
}

} // namespace opensn
