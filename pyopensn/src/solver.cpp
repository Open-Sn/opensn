// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "pyapi.hpp"

#include <cstddef>
#include <memory>

#include "framework/field_functions/field_function_grid_based.h"
#include "framework/physics/solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"

namespace opensn {

// Wrap solver
static void wrap_solver(py::module & slv) {
    // multi-group cross section
    auto solver = py::class_<Solver, std::shared_ptr<Solver>>(
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
        [](LBSSolver & self) {
            py::list field_function_list_per_group;
            for (std::size_t group = 0; group < self.GetNumGroups(); group++) {
                py::list field_function_list_per_moment;
                for (std::size_t moment = 0; moment < self.GetNumMoments(); moment++) {
                    std::size_t ff_index = self.MapPhiFieldFunction(group, moment);
                    field_function_list_per_moment.append(self.GetFieldFunctions()[ff_index]);
                }
                field_function_list_per_group.append(field_function_list_per_moment);
            }
            return field_function_list_per_group;
        },
        R"(
        Return, for each group, a list of field functions corresponding to each moment. Note that the moment index
        varies more rapidly than the group index.
        )"
    );
    // discrete ordinate solver
    auto discrete_ordinate_solver = py::class_<DiscreteOrdinatesSolver, std::shared_ptr<DiscreteOrdinatesSolver>, LBSSolver>(
        slv,
        "DiscreteOrdinatesSolver",
        R"(
        Base class for discrete ordinates solvers.

        This class mostly establishes utilities related to sweeping. From here, steady-state, transient, adjoint, and
        k-eigenvalue solver can be derived.

        Wrapper of :cpp:class:`opensn::DiscreteOrdinatesSolver`.
        )"
    );
    discrete_ordinate_solver.def(
        py::init(
            [](py::kwargs & params) {
                return DiscreteOrdinatesSolver::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a discrete ordinates solver object.
        )"
    );
    // steady dtate solver
    auto steady_state_solver = py::class_<SteadyStateSolver, std::shared_ptr<SteadyStateSolver>, Solver>(
        slv,
        "SteadyStateSolver",
        R"(
        ...

        Wrapper of :cpp:class:`opensn::SteadyStateSolver`.
        )"
    );
    steady_state_solver.def(
        py::init(
            [](py::kwargs & params) {
                return SteadyStateSolver::Create(kwargs_to_param_block(params));
            }
        ),
        R"(
        Construct a steady state solver object.
        )"
    );
}

// Wrap the solver components of OpenSn
void py_solver(py::module & pyopensn) {
    py::module slv = pyopensn.def_submodule("solver", "Solver module.");
    wrap_solver(slv);
}

}  // namespace opensn
