#include "diffusion_solver.h"

#include "chi_log.h"

/**\defgroup LuaDiffusionBasicOptions Basic Options
 * \ingroup LuaDiffusion
 *
Option name           | Type   | Default Value | Description
----------------------|------- |---------------|------------
"discretization_method" | string | "None"        | Spatial discretization method. "PWLC", "PWLD_MIP".
"max_iters"             | int    | 500           | Maximum iterations for solver convergence.
"residual_tolerance"    | float  | 1.0e-8        | Residual convergence tolerance.
"property_map_D"        | int    | 0             | Material property index to use for diffusion coefficient
"property_map_q"        | int    | 1             | Material property index to use for source.
"property_map_sigma"    | int    | 2             | Material property index to use for interaction coefficient.
 *
 * To set these options use the command chiSolverSetBasicOption() with
 * an option-name and value in the table above.
 *
 * ## Example
 * Example usage
 * \code
 * chiSolverSetBasicOption(phys1, "discretization_method", "PWLC")
 * \endcode*/

chi_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, {{"discretization_method", std::string("None")},
                                       {"max_iters", int64_t(500)},
                                       {"residual_tolerance", 1.0e-8},
                                       {"property_map_D",int64_t(0)},
                                       {"property_map_q",int64_t(1)},
                                       {"property_map_sigma",int64_t(2)}})
{}

chi_diffusion::Solver::~Solver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
  KSPDestroy(&ksp_);
}