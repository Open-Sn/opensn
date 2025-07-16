// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/event_system/physics_event_publisher.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/physics/solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/nl_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_scdsa_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_smm_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
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
    [](Solver & self)
    {
      PhysicsEventPublisher::GetInstance().SolverInitialize(self);
    },
    "Initialize the solver."
  );
  solver.def(
    "Execute",
    [](Solver & self)
    {
      PhysicsEventPublisher::GetInstance().SolverExecute(self);
    },
    "Execute the solver."
  );
  solver.def(
    "Step",
    [](Solver & self)
    {
      PhysicsEventPublisher::GetInstance().SolverStep(self);
    },
    "Step the solver."
  );
  solver.def(
    "Advance",
    [](Solver & self)
    {
      PhysicsEventPublisher::GetInstance().SolverAdvance(self);
    },
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
    spatial_discretization: str, default='pwld'
        What spatial discretization to use. Currently only ``pwld`` is supported.
    scattering_order: int, default=1
        The level of harmonic expansion for the scattering source.
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
    adjoint: bool, default=False
        Flag for toggling whether the solver is in adjoint mode.
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
    boundary_conditions: List[Dict], default=[]
        A list containing tables for each boundary specification.
    clear_boundary_conditions: bool, default=False
        Clears all boundary conditions. If no additional boundary conditions are supplied, all
        boundaries become vacuum.
    point_sources: List[pyopensn.source.PointSource], default=[]
        A list of point sources.
    clear_point_sources: bool, default=False
        Clear all point sources.
    volumetric_sources: List[pyopensn.source.VolumetricSource], default=[]
        A list of volumetric sources.
    clear_volumetric_sources: bool, default=False
        Clear all volumetric sources.
    )"
  );
  lbs_problem.def(
    "ComputeFissionRate",
    [](LBSProblem& self, const std::string& scalar_flux_iterate)
    {
      const std::vector<double>* phi_ptr;
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
    [](LBSProblem& self, const std::string& file_base)
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
    [](LBSProblem& self, const std::string& file_base)
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
    "ComputeBalance",
    &DiscreteOrdinatesProblem::ComputeBalance,
    R"(
    Compute and print particle balance for the problem.
    )"
  );
  do_problem.def(
    "ComputeLeakage",
    [](DiscreteOrdinatesProblem& self, py::list bnd_names)
    {
      // get the supported boundaries
      std::map<std::string, std::uint64_t> allowed_bd_names = LBSProblem::supported_boundary_names;
      std::map<std::uint64_t, std::string> allowed_bd_ids = LBSProblem::supported_boundary_ids;
      // get the boundaries to parse
      std::vector<std::uint64_t> bndry_ids;
      if (bnd_names.size() > 1)
      {
        for (py::handle name : bnd_names)
        {
          bndry_ids.push_back(allowed_bd_names.at(name.cast<std::string>()));
        }
      }
      else
      {
        bndry_ids = self.GetGrid()->GetUniqueBoundaryIDs();
      }
      // compute the leakage
      std::map<std::uint64_t, std::vector<double>> leakage = self.ComputeLeakage(bndry_ids);
      // convert result to native Python
      py::dict result;
      for (const auto& [bndry_id, gr_wise_leakage] : leakage)
      {
        py::array_t<double> np_vector = py::array_t<double>(gr_wise_leakage.size());
        py::buffer_info buffer = np_vector.request();
        double* np_vector_data = static_cast<double*>(buffer.ptr);
        std::copy(gr_wise_leakage.begin(), gr_wise_leakage.end(), np_vector_data);
        result[allowed_bd_ids.at(bndry_id).data()] = np_vector;
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
  auto steady_state_solver = py::class_<SteadyStateSolver, std::shared_ptr<SteadyStateSolver>,
                                        Solver>(
    slv,
    "SteadyStateSolver",
    R"(
    Steady state solver.

    Wrapper of :cpp:class:`opensn::SteadyStateSolver`.
    )"
  );
  steady_state_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return SteadyStateSolver::Create(kwargs_to_param_block(params));
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
    lbs_problem: pyopensn.solver.LBSProblem
        Existing LBSProblem instance.
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

  // power iteration k-eigen SCDSA solver
  auto pi_k_eigen_scdsa_solver = py::class_<PowerIterationKEigenSCDSASolver,
                                            std::shared_ptr<PowerIterationKEigenSCDSASolver>,
                                            PowerIterationKEigenSolver>(
    slv,
    "PowerIterationKEigenSCDSASolver",
    R"(
    Power iteration k-eigenvalue solver with SCDSA.

    Wrapper of :cpp:class:`opensn::PowerIterationKEigenSCDSASolver`.
    )"
  );
  pi_k_eigen_scdsa_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PowerIterationKEigenSCDSASolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a power iteration k-eigenvalue solver with SCDSA.

    Parameters
    ----------
    lbs_problem: pyopensn.solver.LBSProblem
        Existing LBSProblem instance.
    max_iters: int, default=1000
        Maximum power iterations allowed.
    k_tol: float, default=1.0e-10
        Tolerance on the k-eigenvalue.
    reset_solution: bool, default=True
        If true, initialize flux moments to 1.0.
    reset_phi0: bool, default=True
        If true, reinitializes scalar fluxes to 1.0.
    accel_pi_max_its : int, default=50
        Maximum number of iterations allowed for the inner power iterations used by the
        acceleration method.
    accel_pi_k_tol : float, default=1.0e-10
        Convergence tolerance on the k-eigenvalue within the inner power iterations of the
        acceleration method.
    accel_pi_verbose : bool, default=False
        If True, enables verbose logging output from the acceleration method's inner solver.
    diff_accel_diffusion_l_abs_tol : float, default=1.0e-10
        Absolute residual tolerance for convergence of the diffusion accelerator.
    diff_accel_diffusion_max_iters : int, default=100
        Maximum number of iterations allowed for the diffusion accelerator solve.
    diff_accel_diffusion_verbose : bool, default=False
        If True, enables verbose logging output from the diffusion accelerator.
    diff_accel_diffusion_petsc_options : str, default="ssss"
        Additional PETSc options passed to the diffusion accelerator linear solver.
    diff_accel_sdm : {'pwld', 'pwlc'}, default='pwld'
        Spatial discretization method to use for the diffusion solver. Valid choices are:
            - 'pwld' : Piecewise Linear Discontinuous
            - 'pwlc' : Piecewise Linear Continuous
    )"
  );

  // power iteration k-eigen SMM solver
  auto pi_k_eigen_smm_solver = py::class_<PowerIterationKEigenSMMSolver,
                                          std::shared_ptr<PowerIterationKEigenSMMSolver>,
                                          PowerIterationKEigenSolver>(
    slv,
    "PowerIterationKEigenSMMSolver",
    R"(
    Power iteration k-eigenvalue solver with SMM acceleration.

    Wrapper of :cpp:class:`opensn::PowerIterationKEigenSMMSolver`.
    )"
  );
  pi_k_eigen_smm_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PowerIterationKEigenSMMSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a power iteration k-eigen SMM solver.

    Parameters
    ----------
    lbs_problem: pyopensn.solver.LBSProblem
        Existing LBSProblem instance.
    max_iters: int, defauilt=1000
        Maximum power iterations allowed.
    k_tol: float, default=1.0e-10
        Tolerance on the k-eigenvalue.
    reset_solution: bool, default=True
        If true, initialize flux moments to 1.0.
    reset_phi0: bool, default=True
        If true, reinitializes scalar fluxes to 1.0.
    accel_pi_max_its: int, default=50
        Maximum number of power iterations allowed for the second-moment method's diffusion solver.
    accel_pi_k_tol: float, default=1.0e-10
        Convergence tolerance for the k-eigenvalue in the second-moment method diffusion solver.
    accel_pi_verbose: bool, default=False
        If True, enables verbose output from the second-moment method diffusion solver.
    diff_sdm: {'pwlc', 'pwld'}, default='pwlc'
        Spatial discretization method for the second-moment method diffusion system.
            - 'pwlc': Piecewise Linear Continuous
            - 'pwld': Piecewise Linear Discontinuous
    diff_l_abs_tol: float, default=1.0e-10
        Absolute residual tolerance for convergence of the diffusion accelerator linear solver.
    diff_l_max_its: int, default=100
        Maximum number of iterations allowed for the diffusion accelerator linear solver.
    diff_petsc_options: str, default=""
        Additional PETSc options to pass to the diffusion accelerator solver.
    diff_verbose: bool, default=False
        If True, enables verbose output from the diffusion accelerator.
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
  WrapPIteration(slv);
}

} // namespace opensn
