#include "lbsadj_solver.h"

#include "ChiObjectFactory.h"

namespace lbs
{

RegisterChiObject(lbs, DiscreteOrdinatesAdjointSolver);

// ##################################################################
/**Returns the input parameters.*/
chi::InputParameters
DiscreteOrdinatesAdjointSolver::GetInputParameters()
{
  chi::InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();

  params.SetGeneralDescription("Adjoint capability");

  params.SetClassName("DiscreteOrdinatesAdjointSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesAdjointSolver");

  return params;
}

// ###################################################################
/**Constructor.*/
DiscreteOrdinatesAdjointSolver::DiscreteOrdinatesAdjointSolver(const chi::InputParameters& params)
  : lbs::DiscreteOrdinatesSolver(params)
{
  basic_options_.AddOption<std::string>("REFERENCE_RF", std::string());
}

// ###################################################################
/**Constructor.*/
lbs::DiscreteOrdinatesAdjointSolver::DiscreteOrdinatesAdjointSolver(const std::string& solver_name)
  : lbs::DiscreteOrdinatesSolver(solver_name)
{
  basic_options_.AddOption<std::string>("REFERENCE_RF", std::string());
}

/**Returns the list of volumetric response functions.*/
const std::vector<lbs::DiscreteOrdinatesAdjointSolver::RespFuncAndSubs>&
lbs::DiscreteOrdinatesAdjointSolver::GetResponseFunctions() const
{
  return response_functions_;
}

} // namespace lbs