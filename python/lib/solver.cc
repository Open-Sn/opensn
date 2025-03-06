// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/physics/solver.h"
#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/lbs_mip_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_solver/lbs_curvilinear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"
#include "modules/linear_boltzmann_solvers/executors/nl_keigen.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_scdsa.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_smm.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/io/lbs_solver_io.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
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
void WrapSolver(py::module& slv)
{
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
}

// Wrap LBS solver
void WrapLBS(py::module& slv)
{
  // LBS solver
  auto lbs_solver = py::class_<LBSSolver, std::shared_ptr<LBSSolver>, Solver>(
    slv,
    "LBSSolver",
    R"(
    Base class for all Linear Boltzmann solvers.

    Wrapper of :cpp:class:`opensn::LBSSolver`.
    )"
  );
  lbs_solver.def(
    "GetScalarFieldFunctionList",
    [](LBSSolver& self)
    {
      py::list field_function_list_per_group;
      for (std::size_t group = 0; group < self.GetNumGroups(); group++)
      {
        py::list field_function_list_per_moment;
        for (std::size_t moment = 0; moment < self.GetNumMoments(); moment++)
        {
          std::size_t ff_index = self.MapPhiFieldFunction(group, moment);
          field_function_list_per_moment.append(self.GetFieldFunctions()[ff_index]);
        }
        field_function_list_per_group.append(field_function_list_per_moment);
      }
      return field_function_list_per_group;
    },
    R"(
    Return, for each group, a list of field functions corresponding to each moment. Note that the
    moment index varies more rapidly than the group index.
    )"
  );
  lbs_solver.def(
    "ComputeFissionRate",
    [](LBSSolver& self, const std::string& nature)
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
  lbs_solver.def(
    "WriteFluxMoments",
    [](LBSSolver& self, const std::string& file_base)
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
  lbs_solver.def(
    "CreateAndWriteSourceMoments",
    [](LBSSolver& self, const std::string& file_base)
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
  lbs_solver.def(
    "ReadFluxMomentsAndMakeSourceMoments",
    [](LBSSolver& self, const std::string& file_base, bool single_file_flag)
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
  lbs_solver.def(
    "ReadSourceMoments",
    [](LBSSolver& self, const std::string& file_base, bool single_file_flag)
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
  lbs_solver.def(
    "ReadFluxMoments",
    [](LBSSolver& self, const std::string& file_base, bool single_file_flag)
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
  lbs_solver.def(
    "WriteAngularFluxes",
    [](LBSSolver& self, const std::string& file_base)
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
  lbs_solver.def(
    "ReadAngularFluxes",
    [](LBSSolver& self, const std::string& file_base)
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
  auto do_solver = py::class_<DiscreteOrdinatesSolver, std::shared_ptr<DiscreteOrdinatesSolver>,
                              LBSSolver>(
    slv,
    "DiscreteOrdinatesSolver",
    R"(
    Base class for discrete ordinates solvers.

    This class mostly establishes utilities related to sweeping. From here, steady-state, transient,
    adjoint, and k-eigenvalue solver can be derived.

    Wrapper of :cpp:class:`opensn::DiscreteOrdinatesSolver`.
    )"
  );
  do_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return DiscreteOrdinatesSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a discrete ordinates solver object.

    Parameters
    ----------
    ???
    )"
  );
  do_solver.def(
    "ComputeBalance",
    &DiscreteOrdinatesSolver::ComputeBalance,
    R"(
    ???
    )"
  );
  do_solver.def(
    "ComputeLeakage",
    [](DiscreteOrdinatesSolver& self, py::list bnd_names)
    {
      // get the supported boundaries
      std::map<std::string, std::uint64_t> allowed_bd_names = LBSSolver::supported_boundary_names;
      std::map<std::uint64_t, std::string> allowed_bd_ids = LBSSolver::supported_boundary_ids;
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
    }
  );

  // discrete ordinate curvilinear solver
  auto do_curvilinear_solver = py::class_<DiscreteOrdinatesCurvilinearSolver,
                                          std::shared_ptr<DiscreteOrdinatesCurvilinearSolver>,
                                          DiscreteOrdinatesSolver>(
    slv,
    "DiscreteOrdinatesCurvilinearSolver",
    R"(
    A neutral particle transport solver in point-symmetric and axial-symmetric curvilinear
    coordinates.

    Wrapper of :cpp:class:`opensn::DiscreteOrdinatesCurvilinearSolver`.
    )"
  );
  do_curvilinear_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return DiscreteOrdinatesCurvilinearSolver::Create(kwargs_to_param_block(params));
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
                                          LBSSolver>(
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
}

// Wrap steady-state solver
void WrapSteadyState(py::module& slv)
{
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
}

// Wrap non-linear k-eigen solver
void WrapNLKEigen(py::module& slv)
{
  // non-linear k-eigen solver
  auto non_linear_k_eigen_solver = py::class_<NonLinearKEigen, std::shared_ptr<NonLinearKEigen>,
                                              Solver>(
    slv,
    "NonLinearKEigen",
    R"(
    ???

    Wrapper of :cpp:class:`opensn::NonLinearKEigen`.
    )"
  );
  non_linear_k_eigen_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return NonLinearKEigen::Create(kwargs_to_param_block(params));
      }
        ),
    R"(
    Construct a non-linear k-eigen solver object.

    Parameters
    ----------
    ???
    )"
  );
}

// Wrap power iteration solvers
void WrapPIteration(py::module& slv)
{
  // power iteration k-eigen solver
  auto pi_k_eigen_solver = py::class_<PowerIterationKEigen, std::shared_ptr<PowerIterationKEigen>,
                                      Solver>(
    slv,
    "PowerIterationKEigen",
    R"(
    ???

    Wrapper of :cpp:class:`opensn::PowerIterationKEigen`.
    )"
  );
  pi_k_eigen_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PowerIterationKEigen::Create(kwargs_to_param_block(params));
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
  auto pi_k_eigen_scdsa_solver = py::class_<PowerIterationKEigenSCDSA,
                                            std::shared_ptr<PowerIterationKEigenSCDSA>,
                                            PowerIterationKEigen>(
    slv,
    "PowerIterationKEigenSCDSA",
    R"(
    ???

    Wrapper of :cpp:class:`opensn::PowerIterationKEigenSCDSA`.
    )"
  );
  pi_k_eigen_scdsa_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PowerIterationKEigenSCDSA::Create(kwargs_to_param_block(params));
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
  auto pi_k_eigen_smm_solver = py::class_<PowerIterationKEigenSMM,
                                          std::shared_ptr<PowerIterationKEigenSMM>,
                                          PowerIterationKEigen>(
    slv,
    "PowerIterationKEigenSMM",
    R"(
    ???

    Wrapper of :cpp:class:`opensn::PowerIterationKEigenSMM`.
    )"
  );
  pi_k_eigen_smm_solver.def(
    py::init(
      [](py::kwargs& params)
      {
        return PowerIterationKEigenSMM::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a power iteration k-eigen SMM solver object.

    Parameters
    ----------
    ???
    )"
  );
}

// Wrap PRK solver
void WrapPRK(py::module& slv)
{
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
    )"
  );
}

// Wrap the solver components of OpenSn
void py_solver(py::module& pyopensn)
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
