#include "lbs_discrete_ordinates_solver.h"

#include "ChiObjectFactory.h"

namespace lbs
{
RegisterChiObject(lbs, DiscreteOrdinatesSolver);
}

lbs::DiscreteOrdinatesSolver::DiscreteOrdinatesSolver(
  const std::string& text_name)
  : LBSSolver(text_name)
{
}

chi::InputParameters lbs::DiscreteOrdinatesSolver::GetInputParameters()
{
  chi::InputParameters params = LBSSolver::GetInputParameters();

  params.SetClassName("DiscreteOrdinatesSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "LBSDiscreteOrdinatesSolver");

  params.AddOptionalParameterArray(
    "directions_sweep_order_to_print",
    std::vector<size_t>(),
    "List of direction id's for which sweep ordering info is to be printed.");

  params.AddOptionalParameter(
    "sweep_type", "AAH", "The sweep type to use for sweep operatorations.");

  using namespace chi_data_types;
  params.ConstrainParameterRange("sweep_type",
                                 AllowableRangeList::New({"AAH", "CBC"}));

  return params;
}

/**Static registration based constructor.*/
lbs::DiscreteOrdinatesSolver::DiscreteOrdinatesSolver(
  const chi::InputParameters& params)
  : LBSSolver(params),
    verbose_sweep_angles_(
      params.GetParamVectorValue<size_t>("directions_sweep_order_to_print")),
    sweep_type_(params.GetParamValue<std::string>("sweep_type"))
{
}

/**Destructor for LBS*/
lbs::DiscreteOrdinatesSolver::~DiscreteOrdinatesSolver()
{
  for (auto& groupset : groupsets_)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);
  }
}

/**Gets the local and global number of iterative unknowns. This normally is
 * only the flux moments, however, the sweep based solvers might include
 * delayed angular fluxes in this number.*/
std::pair<size_t, size_t>
lbs::DiscreteOrdinatesSolver::GetNumPhiIterativeUnknowns()
{
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  size_t num_local_psi_dofs = 0;
  size_t num_globl_psi_dofs = 0;
  for (auto& groupset : groupsets_)
  {
    const auto num_delayed_psi_info =
      groupset.angle_agg_->GetNumDelayedAngularDOFs();
    num_local_psi_dofs += num_delayed_psi_info.first;
    num_globl_psi_dofs += num_delayed_psi_info.second;
  }

  const size_t num_local_dofs = num_local_phi_dofs + num_local_psi_dofs;
  const size_t num_globl_dofs = num_globl_phi_dofs + num_globl_psi_dofs;

  return {num_local_dofs, num_globl_dofs};
}