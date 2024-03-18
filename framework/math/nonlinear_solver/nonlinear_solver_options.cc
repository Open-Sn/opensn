#include "framework/math/nonlinear_solver/nonlinear_solver_options.h"

#include "framework/object_factory.h"

namespace opensn
{

InputParameters
NonLinearSolverOptions::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("Options available on NonLinearSolver");
  params.SetDocGroup("LuaMath");

  params.AddOptionalParameter("name", "NonLinearSolver", "A name to assign to the solver.");

  params.AddOptionalParameter("nl_method", "JFNK", "The non-linear method to use.");
  params.AddOptionalParameter("l_method", "gmres", "The linear solver method to use.");

  params.AddOptionalParameter("petsc_snes_type",
                              "newtonls",
                              "The type passed to SNESSetType, if PETSc is used. Consult PETSc's "
                              "documentation for a list of available types.");

  ParameterBlock pc_options;
  pc_options.AddParameter("pc_type", "hypre");
  pc_options.AddParameter("pc_hypre_type", "boomeramg");

  params.AddOptionalParameterBlock(
    "pc_options", pc_options, "A table of parameters used in the preconditioner.");

  params.AddOptionalParameter("nl_rel_tol", 1.0e-8, "Non-linear relative tolerance");
  params.AddOptionalParameter("nl_abs_tol", 1.0e-8, "Non-linear absolute tolerance");
  params.AddOptionalParameter("nl_sol_tol", 1.0e-50, "Non-linear solution tolerance");
  params.AddOptionalParameter("nl_max_its", 50, "Non-linear maximum iterations");
  params.AddOptionalParameter(
    "nl_max_r_evaluations",
    -1,
    "The maximum allowed residual evaluations. Negative number disables this.");
  params.AddOptionalParameter("l_max_failed_iterations",
                              1000,
                              "The maximum allowed non-linear iterations "
                              "where the linear solver failed to converge.");

  params.AddOptionalParameter("l_rel_tol", 1.0e-8, "Linear relative tolerance");
  params.AddOptionalParameter("l_abs_tol", 1.0e-8, "Linear absolute tolerance");
  params.AddOptionalParameter("l_div_tol", 1.0e6, "Linear divergence tolerance");
  params.AddOptionalParameter("l_max_its", 100, "Linear maximum iterations");
  params.AddOptionalParameter("l_gmres_restart_intvl", 30, "GMRes restart interval");
  params.AddOptionalParameter("l_gmres_breakdown_tol", 1.0e6, "GMRes breakdown tolerance");

  params.ConstrainParameterRange("nl_method",
                                 AllowableRangeList::New({"JFNK", "PJFNK", "NEWTON", "LINEAR"}));

  return params;
}

NonLinearSolverOptions::NonLinearSolverOptions(const InputParameters& params)
  : Object(params),
    nl_method_(params.GetParamValue<std::string>("nl_method")),
    l_method_(params.GetParamValue<std::string>("l_method")),
    pc_options_(params.GetParam("pc_options")),
    petsc_snes_type_(params.GetParamValue<std::string>("petsc_snes_type")),
    nl_rel_tol_(params.GetParamValue<double>("nl_rel_tol")),
    nl_abs_tol_(params.GetParamValue<double>("nl_abs_tol")),
    nl_sol_tol_(params.GetParamValue<double>("nl_sol_tol")),
    nl_max_its_(params.GetParamValue<int>("nl_max_its")),
    nl_max_r_evaluations_(params.GetParamValue<int>("nl_max_r_evaluations")),
    l_max_failed_iterations_(params.GetParamValue<int>("l_max_failed_iterations")),
    l_rel_tol_(params.GetParamValue<double>("l_rel_tol")),
    l_abs_tol_(params.GetParamValue<double>("l_abs_tol")),
    l_div_tol_(params.GetParamValue<double>("l_div_tol")),
    l_max_its_(params.GetParamValue<int>("l_max_its")),
    l_gmres_restart_intvl_(params.GetParamValue<int>("l_gmres_restart_intvl")),
    l_gmres_breakdown_tol_(params.GetParamValue<double>("l_gmres_breakdown_tol"))
{
}

} // namespace opensn
