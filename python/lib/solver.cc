// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/event_system/physics_event_publisher.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/physics/solver.h"
#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/lbs_mip_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/lbs_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"
#include "modules/linear_boltzmann_solvers/executors/nl_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_scdsa_solver.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_smm_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/point_reactor_kinetics/point_reactor_kinetics.h"
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
  solver.def(
    "GetFieldFunctions",
    [](Solver& self)
    {
      py::list ff_list;
      std::vector<std::shared_ptr<FieldFunctionGridBased>>& cpp_ff_list = self.GetFieldFunctions();
      for (std::shared_ptr<FieldFunctionGridBased>& ff : cpp_ff_list) {
        ff_list.append(ff);
      }
      return ff_list;
    },
    "Get the list of field functions."
  );
  // clang-format on
}

// Wrap LBS solver
void
WrapLBS(py::module& slv)
{
  // clang-format off
  // LBS problem
  auto lbs_problem = py::class_<LBSProblem, std::shared_ptr<LBSProblem>, Solver>(
    slv,
    "LBSProblem",
    R"(
    Base class for all Linear Boltzmann solvers.

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
    Return, for each group, a list of field functions corresponding to each moment. Note that the
    moment index varies more rapidly than the group index.

    Parameters
    ----------
    only_scalar_flux: bool, default=True
        If True, only return a list of of field functions corresponding to moment zero-th for each
        group. The result is only a simple list of field functions. Otherwise, the result will be a
        list per group of list per moment.
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
    Set options from a large list of parameters.

    Parameters
    ----------
    spatial_discretization: str, default='pwld'
        What spatial discretization to use. Currently only ``pwld`` is supported.
    scattering_order: int, default=1
        The level of harmonic expansion for the scattering source.
    max_mpi_message_size: int default=32768
        The maximum MPI message size used during sweep initialization.
    restart_writes_enabled: bool, default=False
        Flag that controls writing of restart dumps.
    write_delayed_psi_to_restart: bool, default=True
        Flag that controls writing of delayed angular fluxes to restarts.
    read_restart_path: str, default=''
        Full path for reading restart dumps including file stem.
    write_restart_path: str, default=''
        Full path for writing restart dumps including file stem.
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
    [](LBSProblem& self, const std::string& nature)
    {
      const std::vector<double>* phi_ptr;
      if (nature == "old")
      {
        phi_ptr = &self.GetPhiOldLocal();
      }
      else if (nature == "new")
      {
        phi_ptr = &self.GetPhiNewLocal();
      }
      else
      {
        throw std::invalid_argument("Unknown nature type \"" + nature + "\".");
      }
      return self.ComputeFissionRate(*phi_ptr);
    },
    R"(
    ???

    Parameters
    ----------
    nature: {'old', 'new'}
        ???
    )",
    py::arg("nature")
  );
  lbs_problem.def(
    "WriteFluxMoments",
    [](LBSProblem& self, const std::string& file_base)
    {
      LBSSolverIO::WriteFluxMoments(self, file_base);
    },
    R"(
    ???

    Parameters
    ----------
    file_base: str
        ???
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
    ???

    Parameters
    ----------
    file_base: str
        ???
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
    ???

    Parameters
    ----------
    file_base: str
        ???
    single_file_flag: bool
        ???
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
    ???

    Parameters
    ----------
    file_base: str
        ???
    single_file_flag: bool
        ???
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
    ???

    Parameters
    ----------
    file_base: str
        ???
    single_file_flag: bool
        ???
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
    ???

    Parameters
    ----------
    file_base: str
        ???
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
    ???

    Parameters
    ----------
    file_base: str
        ???
    )",
    py::arg("file_base")
  );

  // discrete ordinate solver
  auto do_problem = py::class_<DiscreteOrdinatesProblem, std::shared_ptr<DiscreteOrdinatesProblem>,
                               LBSProblem>(
    slv,
    "DiscreteOrdinatesProblem",
    R"(
    Base class for discrete ordinates solvers.

    This class mostly establishes utilities related to sweeping. From here, steady-state, transient,
    adjoint, and k-eigenvalue solver can be derived.

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
    Construct a discrete ordinates solver object.

    Parameters
    ----------
    ???
    )"
  );
  do_problem.def(
    "ComputeBalance",
    &DiscreteOrdinatesProblem::ComputeBalance,
    R"(
    ???
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
    ???

    Parameters
    ----------
    bnd_names: List[str]
        ???
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
    A neutral particle transport solver in point-symmetric and axial-symmetric curvilinear
    coordinates.

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
    Construct a discrete curvilinear ordinates solver object.

    Parameters
    ----------
    ???
    )"
  );

  // diffusion DFEM solver
  auto diffusion_dfem_solver = py::class_<DiffusionDFEMSolver, std::shared_ptr<DiffusionDFEMSolver>,
                                          LBSProblem>(
    slv,
    "DiffusionDFEMSolver",
    R"(
    ???

    Wrapper of :cpp:class:`opensn::DiffusionDFEMSolver`.
    )"
  );
  diffusion_dfem_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return DiffusionDFEMSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a diffusion DFEM solver object.

    Parameters
    ----------
    ???
    )"
  );
  // clang-format on
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
    ???

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
    Construct a steady state solver object.

    Parameters
    ----------
    ???
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
    ???

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
    Construct a non-linear k-eigen solver object.

    Parameters
    ----------
    ???
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
    ???

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
    Construct a power iteration k-eigen solver object.

    Parameters
    ----------
    ???
    )"
  );

  // power iteration k-eigen SCDSA solver
  auto pi_k_eigen_scdsa_solver = py::class_<PowerIterationKEigenSCDSASolver,
                                            std::shared_ptr<PowerIterationKEigenSCDSASolver>,
                                            PowerIterationKEigenSolver>(
    slv,
    "PowerIterationKEigenSCDSASolver",
    R"(
    ???

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
    Construct a power iteration k-eigen SCDSA solver object.

    Parameters
    ----------
    ???
    )"
  );

  // power iteration k-eigen SMM solver
  auto pi_k_eigen_smm_solver = py::class_<PowerIterationKEigenSMMSolver,
                                          std::shared_ptr<PowerIterationKEigenSMMSolver>,
                                          PowerIterationKEigenSolver>(
    slv,
    "PowerIterationKEigenSMMSolver",
    R"(
    ???

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
    Construct a power iteration k-eigen SMM solver object.

    Parameters
    ----------
    ???
    )"
  );
  // clang-format on
}

// Wrap PRK solver
void
WrapPRK(py::module& slv)
{
  // clang-format off
  // point reactor kineatic solver
  auto prk_solver = py::class_<PRKSolver, std::shared_ptr<PRKSolver>, Solver>(
    slv,
    "PRKSolver",
    R"(
    General transient solver for point kinetics.

    Wrapper of :cpp:class:`opensn::PRKSolver`.
    )"
  );
  prk_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PRKSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a point reactor kineatic solver object.

    Parameters
    ----------
    ???
    )"
  );
  prk_solver.def(
    "GetPopulationPrev",
    &PRKSolver::GetPopulationPrev,
    R"(
    Get the population at the previous time step.
    )"
  );
  prk_solver.def(
    "GetPopulationNew",
    &PRKSolver::GetPopulationNew,
    R"(
    Get the population at the next time step.
    )"
  );
  prk_solver.def(
    "GetPeriod",
    &PRKSolver::GetPeriod,
    R"(
    Get the period computed for the last time step.
    )"
  );
  prk_solver.def(
    "GetTimePrev",
    &PRKSolver::GetTimePrev,
    R"(
    Get the time computed for the last time step.
    )"
  );
  prk_solver.def(
    "GetTimeNew",
    &PRKSolver::GetTimeNew,
    R"(
    Get the time computed for the next time step.
    )"
  );
  prk_solver.def(
    "GetSolutionPrev",
    [](PRKSolver& self)
    {
      return convert_vector(self.GetSolutionPrev());
    },
    R"(
    Get the solution at the previous time step.
    )"
  );
  prk_solver.def(
    "GetSolutionNew",
    [](PRKSolver& self)
    {
      return convert_vector(self.GetSolutionNew());
    },
    R"(
    Get the solution at the next time step.
    )"
  );
  prk_solver.def(
    "SetRho",
    &PRKSolver::SetRho,
    R"(
    Set the value of rho.
    ??? (what is rho?)
    )",
    py::arg("rho")
  );
  // clang-format on
}

// Wrap the solver components of OpenSn
void
py_solver(py::module& pyopensn)
{
  py::module slv = pyopensn.def_submodule("solver", "Solver module.");
  WrapSolver(slv);
  WrapLBS(slv);
  WrapSteadyState(slv);
  WrapNLKEigen(slv);
  WrapPIteration(slv);
  WrapPRK(slv);
}

} // namespace opensn
