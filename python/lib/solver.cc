// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/runtime.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/discrete_ordinates_keigen_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/scdsa_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/smm_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/nl_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
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
  problem.def(
    "GetFieldFunctions",
    [](Problem& self)
    {
      py::list ff_list;
      std::vector<std::shared_ptr<FieldFunctionGridBased>>& cpp_ff_list = self.GetFieldFunctions();
      for (std::shared_ptr<FieldFunctionGridBased>& ff : cpp_ff_list) {
        ff_list.append(ff);
      }
      return ff_list;
    },
    R"(
    Get the list of field functions.

    Returns
    -------
    List[pyopensn.fieldfunc.FieldFunctionGridBased]
        List of grid-based field functions representing solution data such as scalar fluxes.
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
    "Step",
    &Solver::Step,
    "Step the solver."
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
    "GetScalarFieldFunctionList",
    [](LBSProblem& self, bool only_scalar_flux)
    {
      py::list field_function_list_per_group;
      for (std::size_t group = 0; group < self.GetNumGroups(); group++)
      {
        if (only_scalar_flux)
        {
          std::size_t ff_index = self.MapPhiFieldFunction(group, 0);
          field_function_list_per_group.append(self.GetFieldFunctions()[ff_index]);
        }
        else
        {
          py::list field_function_list_per_moment;
          for (std::size_t moment = 0; moment < self.GetNumMoments(); moment++)
          {
            std::size_t ff_index = self.MapPhiFieldFunction(group, moment);
            field_function_list_per_moment.append(self.GetFieldFunctions()[ff_index]);
          }
          field_function_list_per_group.append(field_function_list_per_moment);
        }
      }
      return field_function_list_per_group;
    },
    R"(
    Return field functions grouped by energy group and, optionally, by moment.

    Parameters
    ----------
    only_scalar_flux : bool, default=True
        If True, returns only the zeroth moment (scalar flux) field function for each group.
        The result is a flat list of field functions, one per group.

        If False, returns all moment field functions for each group.
        The result is a nested list where each entry corresponds to a group and contains
        a list of field functions for all moments (e.g., scalar flux, higher-order moments).

    Returns
    -------
    Union[List[pyopensn.fieldfunc.FieldFunctionGridBased], List[List[pyopensn.fieldfunc.FieldFunctionGridBased]]]
        The structure of the returned list depends on the `only_scalar_flux` flag.

    Notes
    -----
    The moment index varies more rapidly than the group index when `only_scalar_flux` is False.
    )",
    py::arg("only_scalar_flux") = true
  );
  lbs_problem.def(
    "GetPowerFieldFunction",
    &LBSProblem::GetPowerFieldFunction,
    R"(
    Returns the power generation field function, if enabled.
    )"
  );
  lbs_problem.def(
    "SetOptions",
    [](LBSProblem& self, py::kwargs& params)
    {
      InputParameters input = LBSProblem::GetOptionsBlock();
      input.AssignParameters(kwargs_to_param_block(params));
      self.SetOptions(input);
    },
    R"(
    Set problem options from a large list of parameters.

    Parameters
    ----------
    max_mpi_message_size: int default=32768
        The maximum MPI message size used during sweeps.
    restart_writes_enabled: bool, default=False
        Flag that controls writing of restart dumps.
    write_delayed_psi_to_restart: bool, default=True
        Flag that controls writing of delayed angular fluxes to restarts.
    read_restart_path: str, default=''
        Full path for reading restart dumps including file basename.
    write_restart_path: str, default=''
        Full path for writing restart dumps including file basename.
    write_restart_time_interval: int, default=0
        Time interval in seconds at which restart data is to be written.
    use_precursors: bool, default=False
        Flag for using delayed neutron precursors.
    use_source_moments: bool, default=False
        Flag for ignoring fixed sources and selectively using source moments obtained elsewhere.
    save_angular_flux: bool, default=False
        Flag indicating whether angular fluxes are to be stored or not.
    verbose_inner_iterations: bool, default=True
        Flag to control verbosity of inner iterations.
    verbose_outer_iterations: bool, default=True
        Flag to control verbosity of across-groupset iterations.
    max_ags_iterations: int, default=100
        Maximum number of across-groupset iterations.
    ags_tolerance: float, default=1.0e-6
        Across-groupset iterations tolerance.
    ags_convergence_check: {'l2', 'pointwise'}, default='l2'
        Type of convergence check for AGS iterations.
    verbose_ags_iterations: bool, default=True
        Flag to control verbosity of across-groupset iterations.
    power_field_function_on: bool, default=False
        Flag to control the creation of the power generation field function. If set to ``True``, a
        field function will be created with the general name ``<solver_name>_power_generation``.
    power_default_kappa: float, default=3.20435e-11
        Default ``kappa`` value (Energy released per fission) to use for power generation when cross
        sections do not have ``kappa`` values. Default corresponds to 200 MeV per fission.
    power_normalization: float, default=-1.0
        Power normalization factor to use. Supply a negative or zero number to turn this off.
    field_function_prefix_option: {'prefix', 'solver_name'}, default='prefix'
        Prefix option on field function names. If unset, flux field functions will be exported as
        ``phi_gXXX_mYYY``, where ``XXX`` is the zero-padded 3-digit group number and ``YYY`` is the
        zero-padded 3-digit moment.
    field_function_prefix: str, default=''
        Prefix to use on all field functions. By default, this is empty. If specified, flux moments
        are exported as ``prefix_phi_gXXX_mYYY``.
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
      LBSSolverIO::ReadFluxMoments(self, file_base, single_file_flag, self.GetExtSrcMomentsLocal());
      log.Log() << "Making source moments from flux file.";
      std::vector<double>& temp_phi = self.GetPhiOldLocal();
      self.GetPhiOldLocal() = self.GetExtSrcMomentsLocal();
      self.GetExtSrcMomentsLocal() = self.MakeSourceMomentsFromPhi();
      self.GetPhiOldLocal() = temp_phi;
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
      LBSSolverIO::ReadFluxMoments(self, file_base, single_file_flag, self.GetExtSrcMomentsLocal());
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
      LBSSolverIO::WriteAngularFluxes(self, file_base);
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
      LBSSolverIO::ReadAngularFluxes(self, file_base);
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
    "WriteSurfaceAngularFluxes",
    [](DiscreteOrdinatesProblem& self, 
      const std::string& file_base, 
      py::list bndry_names,
      py::object surfaces)
    {
      // Map boundary names
      std::map<std::string, uint64_t> allowed_bd_names = LBSProblem::supported_boundary_names;
      std::map<std::uint64_t, std::string> allowed_bd_ids = LBSProblem::supported_boundary_ids;
      std::vector<std::string> bndrys;
      for (py::handle name : bndry_names)
        bndrys.push_back(name.cast<std::string>());

      // Map surface names
      std::optional<std::pair<std::string, double>> opt_surfaces;
      if (!surfaces.is_none())
      {
        py::list surf_list = surfaces.cast<py::list>();
        if (py::len(surf_list) == 2)
        {
          std::string surf_id = surf_list[0].cast<std::string>();
          double slice = surf_list[1].cast<double>();
          opt_surfaces = std::make_pair(surf_id, slice);
        }
      }

      // LBSSolverIO::WriteSurfaceAngularFluxes(self, file_base, bndry_map);
      LBSSolverIO::WriteSurfaceAngularFluxes(self, file_base, bndrys, opt_surfaces);
    },
    R"(
    Write surface angular flux data to file.

    Parameters
    ----------
    file_base: str
        File basename.
    )",
    py::arg("file_base"),
    py::arg("bndry_names"),
    py::arg("surfaces") = py::none()
  );
  lbs_problem.def(
    "ReadSurfaceAngularFluxes",
    [](DiscreteOrdinatesProblem& self, const std::string& file_base, py::list bndry_names)
    {
      std::map<std::string, std::uint64_t> supported_bd_names = LBSProblem::supported_boundary_names;
      std::map<std::uint64_t, std::string> supported_bd_ids = LBSProblem::supported_boundary_ids;

      // std::map<std::string, std::uint64_t> bndry_map;
      std::vector<std::string> bndrys;
      for (py::handle name : bndry_names)
      {
        std::string bndry = name.cast<std::string>();
        // bndry_map[bndry] = supported_bd_names.at(bndry);
        bndrys.push_back(bndry);
      }
      // LBSSolverIO::ReadSurfaceAngularFluxes(self, file_base, bndry_map);
      LBSSolverIO::ReadSurfaceAngularFluxes(self, file_base, bndrys);
    },
    R"(
    Read surface angular fluxes from file.

    Parameters
    ----------
    file_base: str
        File basename.
    )",
    py::arg("file_base"),
    py::arg("bndry_names")
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
    "SetBoundaryOptions",
    [](LBSProblem& self, py::kwargs& params)
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
            InputParameters input = LBSProblem::GetBoundaryOptionsBlock();
            input.AssignParameters(pyobj_to_param_block("", boundary.cast<py::dict>()));
            self.SetBoundaryOptions(input);
          }
        }
        else
          throw std::runtime_error("Invalid argument provided to SetBoundaryOptions.\n");
      }
    }
  );
  lbs_problem.def(
    "SetAdjoint",
    [](LBSProblem& self, bool adjoint)
    {
      self.SetAdjoint(adjoint);
    }
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
        groupset.
    xs_map : List[Dict], default=[]
        A list of mappings from block ids to cross-section definitions.
    boundary_conditions: List[Dict], default=[]
        A list containing tables for each boundary specification.
    point_sources: List[pyopensn.source.PointSource], default=[]
        A list of point sources.
    volumetric_sources: List[pyopensn.source.VolumetricSource], default=[]
        A list of volumetric sources.
    options : Dict, default={}
        A block of optional configuration parameters. See `SetOptions` for available settings.
    sweep_type : str, default="AAH"
        The sweep type to use. Must be one of `AAH` or `CBC`. Defaults to `AAH`.
    use_gpus : bool, default=False
        A flag specifying whether GPU acceleration is used for the sweep. Currently, only ``AAH`` is
        supported.
    )"
  );
  do_problem.def(
    "SetOptions",
    [](DiscreteOrdinatesProblem& self, py::kwargs& params)
    {
      InputParameters input = DiscreteOrdinatesProblem::GetOptionsBlock();
      input.AssignParameters(kwargs_to_param_block(params));
      self.SetOptions(input);
    },
    R"(
    Set problem options from a large list of parameters.

    Parameters
    ----------
    adjoint: bool, default=False
        Flag for toggling whether the solver is in adjoint mode.
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
    "ComputeBalance",
    [](DiscreteOrdinatesProblem& self)
    {
      ComputeBalance(self);
    },
    R"(
    Compute and print particle balance for the problem.
    )"
  );
  do_problem.def(
    "ComputeLeakage",
    [](DiscreteOrdinatesProblem& self, py::list bnd_names)
    {
      auto grid = self.GetGrid();
      // get the supported boundaries
      std::map<std::string, std::uint64_t> allowed_bd_names = grid->GetBoundaryNameMap();
      std::map<std::uint64_t, std::string> allowed_bd_ids = grid->GetBoundaryIDMap();
      // get the boundaries to parse, preserving user order
      std::vector<std::uint64_t> bndry_ids;
      if (bnd_names.size() > 1)
      {
        for (py::handle name : bnd_names)
        {
          auto sname = name.cast<std::string>();
          bndry_ids.push_back(allowed_bd_names.at(sname));
        }
      }
      else
      {
        bndry_ids = self.GetGrid()->GetUniqueBoundaryIDs();
      }
      // compute the leakage
      std::map<std::uint64_t, std::vector<double>> leakage = ComputeLeakage(self, bndry_ids);
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
        const std::string& name = allowed_bd_ids.at(bndry_id);
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
        groupset.
    xs_map : list of dict
        A list of mappings from block ids to cross-section definitions.
    boundary_conditions: List[Dict], default=[]
        A list containing tables for each boundary specification.
    point_sources: List[pyopensn.source.PointSource], default=[]
        A list of point sources.
    volumetric_sources: List[pyopensn.source.VolumetricSource], default=[]
        A list of volumetric sources.
    options : dict, optional
        A block of optional configuration parameters. See `SetOptions` for available settings.
    sweep_type : str, optional
        The sweep type to use. Must be one of `AAH` or `CBC`. Defaults to `AAH`.
    )"
  );
}

// Wrap steady-state solver
void
WrapSteadyState(py::module& slv)
{
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
    pyopensn.solver.LBSProblem : LBSProblem
        Existing LBSProblem instance.
    )"
  );
  // clang-format on
}

// Wrap non-linear k-eigen solver
void
WrapNLKEigen(py::module& slv)
{
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
    lbs_problem: pyopensn.solver.LBSProblem
        Existing LBSProblem instance.
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
    Return the current k‑eigenvalue.
    )"
  );
  // clang-format on
}

// Wrap power iteration solvers
void
WrapPIteration(py::module& slv)
{
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
    problem: pyopensn.solver.LBSProblem
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
    )"
  );
  pi_k_eigen_solver.def(
    "GetEigenvalue",
    &PowerIterationKEigenSolver::GetEigenvalue,
    R"(
    Return the current k‑eigenvalue.
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
    problem: pyopensn.solver.LBSProblem
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
    problem: pyopensn.solver.LBSProblem
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
  WrapNLKEigen(slv);
  WrapDiscreteOrdinatesKEigenAcceleration(slv);
  WrapPIteration(slv);
}

} // namespace opensn
