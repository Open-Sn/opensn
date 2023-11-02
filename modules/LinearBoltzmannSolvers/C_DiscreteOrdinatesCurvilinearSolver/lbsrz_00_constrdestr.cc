#include "lbs_curvilinear_solver.h"

#include "ChiObjectFactory.h"

namespace lbs
{

RegisterChiObject(lbs, DiscreteOrdinatesCurvilinearSolver);

chi::InputParameters DiscreteOrdinatesCurvilinearSolver::GetInputParameters()
{
  chi::InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();

  params.SetGeneralDescription(
    "Solver for Discrete Ordinates in cylindrical and spherical coordinates");

  params.SetClassName("DiscreteOrdinatesCurvilinearSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name",
                                       "DiscreteOrdinatesCurvilinearSolver");
  params.AddRequiredParameter<int>("coord_system",
                                   "Coordinate system to use. Can only be 2 or "
                                   "3. 2=Cylindrical, 3=Spherical.");

  return params;
}

DiscreteOrdinatesCurvilinearSolver::DiscreteOrdinatesCurvilinearSolver(
  const chi::InputParameters& params)
  : DiscreteOrdinatesSolver(params),
    coord_system_type_(static_cast<chi_math::CoordinateSystemType>(
      params.GetParamValue<int>("coord_system")))
{
}

} // namespace lbs