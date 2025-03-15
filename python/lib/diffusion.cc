// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "python/lib/functor.h" // temporary, see the included header for more details!
#include "framework/physics/solver.h"
#include "modules/diffusion/cfem_diffusion_solver.h"
#include "modules/diffusion/dfem_diffusion_solver.h"
#include "modules/diffusion/diffusion_solver.h"
#include <memory>

namespace opensn
{

// Wrap diffusion solvers
void
WrapDiffusion(py::module& diffusion)
{
  // clang-format off
  // diffusion solver
  auto diff_base = py::class_<DiffusionSolverBase, std::shared_ptr<DiffusionSolverBase>, Solver>(
    diffusion,
    "DiffusionSolverBase",
    R"(
    Base class for diffusion solvers.

    Wrapper of :cpp:class:`opensn::DiffusionSolverBase`.
    )"
  );
  diff_base.def(
    "UpdateFieldFunctions",
    &DiffusionSolverBase::UpdateFieldFunctions,
    R"(
    Updates the field functions with the latest data.
    )"
  );
  diff_base.def(
    "SetDCoefFunction",
    [](DiffusionSolverBase& self, std::shared_ptr<PySSMFunction> p_func)
    {
      self.SetDCoefFunction(p_func);
    },
    R"(
    ???

    Parameters
    ----------
    func: pyopensn.math.ScalarSpatialMaterialFunction
        ???
    )",
    py::arg("func")
  );
  diff_base.def(
    "SetQExtFunction",
    [](DiffusionSolverBase& self, std::shared_ptr<PySSMFunction> p_func)
    {
      self.SetQExtFunction(p_func);
    },
    R"(
    ???

    Parameters
    ----------
    func: pyopensn.math.ScalarSpatialMaterialFunction
        ???
    )",
    py::arg("func")
  );
  diff_base.def(
    "SetSigmaAFunction",
    [](DiffusionSolverBase& self, std::shared_ptr<PySSMFunction> p_func)
    {
      self.SetSigmaAFunction(p_func);
    },
    R"(
    ???

    Parameters
    ----------
    func: pyopensn.math.ScalarSpatialMaterialFunction
        ???
    )",
    py::arg("func")
  );

  // CFEM diffusion solver
  auto cfem_diffusion = py::class_<CFEMDiffusionSolver, std::shared_ptr<CFEMDiffusionSolver>,
                                   DiffusionSolverBase>(
    diffusion,
    "CFEMDiffusionSolver",
    R"(
    CFEM diffusion solver.

    Diffusion solver using continuous finite-element method.

    Wrapper of :cpp:class:`opensn::CFEMDiffusionSolver`.
    )"
  );
  cfem_diffusion.def(
    py::init(
      [](py::kwargs& params)
      {
        return CFEMDiffusionSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a CFEM diffusion solver.

    Parameters
    ----------
    ???
    )"
  );

  // DFEM diffusion solver
  auto dfem_diffusion = py::class_<DFEMDiffusionSolver, std::shared_ptr<DFEMDiffusionSolver>,
                                   DiffusionSolverBase>(
    diffusion,
    "DFEMDiffusionSolver",
    R"(
    DFEM diffusion solver.

    Diffusion solver using discontinuous finite-element method.

    Wrapper of :cpp:class:`opensn::DFEMDiffusionSolver`.
    )"
  );
  dfem_diffusion.def(
    py::init(
      [](py::kwargs& params)
      {
        return DFEMDiffusionSolver::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a DFEM diffusion solver.

    Parameters
    ----------
    ???
    )"
  );
  // clang-format on
}

// Wrap the diffusion components of OpenSn
void
py_diffusion(py::module& pyopensn)
{
  py::module diffusion = pyopensn.def_submodule("diffusion", "Diffusion function module.");
  WrapDiffusion(diffusion);
}

} // namespace opensn
