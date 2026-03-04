// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#include "framework/object_factory.h"
#include "framework/data_types/vector_ghost_communicator/vector_ghost_communicator.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/runtime.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "modules/diffusion/diffusion_pwlc_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/scdsa_acceleration.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"
namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, SCDSAAcceleration);

InputParameters
SCDSAAcceleration::GetInputParameters()
{
  auto params = DiscreteOrdinatesKEigenAcceleration::GetInputParameters();

  params.AddOptionalParameter("sdm", "pwld", "Spatial discretization");

  params.ConstrainParameterRange("sdm", AllowableRangeList::New({"pwld", "pwlc"}));

  params.ChangeExistingParamToOptional("name", "SCDSAAcceleration");

  return params;
}

std::shared_ptr<SCDSAAcceleration>
SCDSAAcceleration::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SCDSAAcceleration>("lbs::SCDSAAcceleration", params);
}

SCDSAAcceleration::SCDSAAcceleration(const InputParameters& params)
  : DiscreteOrdinatesKEigenAcceleration(params), sdm_(params.GetParamValue<std::string>("sdm"))
{
}

void
SCDSAAcceleration::Initialize()
{
  front_wgs_solver_ = do_problem_.GetWGSSolver(front_gs_.id);
  front_wgs_context_ = std::dynamic_pointer_cast<WGSContext>(front_wgs_solver_->GetContext());
  OpenSnLogicalErrorIf(not front_wgs_context_, ": Casting failed");

  // Make UnknownManager
  const auto num_gs_groups = front_gs_.GetNumGroups();
  UnknownManager uk_man;
  uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

  // Make boundary conditions
  const auto bcs = TranslateBCs(do_problem_.GetSweepBoundaries(), true);

  // Make xs map
  const auto matid_2_mgxs_map =
    PackGroupsetXS(do_problem_.GetBlockID2XSMap(), front_gs_.first_group, front_gs_.last_group);

  // Create solver
  const auto& sdm = do_problem_.GetSpatialDiscretization();
  const auto& unit_cell_matrices = do_problem_.GetUnitCellMatrices();

  if (sdm_ == "pwld")
  {
    diffusion_solver_ = std::make_shared<DiffusionMIPSolver>(std::string(GetName() + "_WGDSA"),
                                                             sdm,
                                                             uk_man,
                                                             bcs,
                                                             matid_2_mgxs_map,
                                                             unit_cell_matrices,
                                                             false,
                                                             true);
  }
  else
  {
    const auto& sdm = do_problem_.GetSpatialDiscretization();
    pwlc_ptr_ = PieceWiseLinearContinuous::New(sdm.GetGrid());
    ghost_info_ = MakePWLDGhostInfo(*pwlc_ptr_, do_problem_.GetUnknownManager());
    diffusion_solver_ = std::make_shared<DiffusionPWLCSolver>(std::string(GetName() + "_WGDSA"),
                                                              *pwlc_ptr_,
                                                              uk_man,
                                                              bcs,
                                                              matid_2_mgxs_map,
                                                              unit_cell_matrices,
                                                              false,
                                                              true);
  }

  auto& ds = *diffusion_solver_;
  ds.options.residual_tolerance = l_abs_tol_;
  ds.options.max_iters = max_iters_;
  ds.options.verbose = verbose_;
  ds.options.additional_options_string = petsc_options_;

  log.Log() << "Initializing diffusion solver";
  ds.Initialize();
  opensn::mpi_comm.barrier();
  log.Log() << "Done Initializing diffusion solver";

  log.Log() << "Assembling A and b";
  std::vector<double> dummy_rhs;
  if (sdm_ == "pwld")
    dummy_rhs.assign(sdm.GetNumLocalDOFs(uk_man), 0.0);
  else
    dummy_rhs.assign(pwlc_ptr_->GetNumLocalAndGhostDOFs(uk_man), 0.0);
  ds.AssembleAand_b(dummy_rhs);
  log.Log() << "Done Assembling A and b";
}

void
SCDSAAcceleration::PreExecute()
{
  phi_temp_ = phi_old_local_;

  // If using multiple sweeps, the algorithm requires a delta-phi between two
  // successive iterations. This is not possible to obtain with the standard solve
  // routine.  Instead, a follow-up sweep must be performed.  If using two or more
  // sweeps, reduce the sweep count by one to ensure the user gets the requested
  // number of sweeps per outer.
  extra_sweep_ = front_gs_.max_iterations > 1;
  if (extra_sweep_)
    front_gs_.max_iterations -= 1;
}

void
SCDSAAcceleration::PrePowerIteration()
{
  Sf_ell_ = q_moments_local_;
  CopyOnlyPhi0(Sf_ell_, Sf0_ell_);

  // Init old scalar flux storage. If using a single sweep, set
  // it to the existing scalar flux
  if (!extra_sweep_)
    CopyOnlyPhi0(phi_old_local_, phi0_ell_);
}

double
SCDSAAcceleration::PostPowerIteration()
{
  // Lambda for the creation of scattering sources but the
  // input vector is only the zeroth moment
  auto SetLBSScatterSourcePhi0 = [this](const std::vector<double>& input,
                                        const bool additive,
                                        const bool suppress_wg_scat = false)
  {
    ProjectBackPhi0(input, phi_temp_);
    solver_->SetLBSScatterSource(phi_temp_, additive, suppress_wg_scat);
  };

  if (extra_sweep_)
  {
    // Set the old scalar flux to the scalar flux before the last sweep
    CopyOnlyPhi0(phi_new_local_, phi0_ell_);

    // Do an extra sweep
    q_moments_local_ = Sf_ell_;
    solver_->SetLBSScatterSource(phi_old_local_, true);
    front_wgs_context_->ApplyInverseTransportOperator(SourceFlags());
  }

  // Store the intermediate scalar flux
  CopyOnlyPhi0(phi_new_local_, phi0_star_);

  // Power Iteration Acceleration
  SetLBSScatterSourcePhi0(phi0_star_ - phi0_ell_, false);
  CopyOnlyPhi0(q_moments_local_, Ss0_res_);

  double production_k = ComputeFissionProduction(do_problem_, phi_new_local_);

  for (const auto& vec : {&epsilon_k_, &epsilon_kp1_})
  {
    vec->resize(phi0_star_.size());
    std::fill(vec->begin(), vec->end(), 0.0);
  }

  double lambda_k = solver_->GetEigenvalue();
  double lambda_kp1 = lambda_k;

  for (int k = 0; k < pi_max_its_; ++k)
  {
    ProjectBackPhi0(epsilon_k_ + phi0_star_, phi_temp_);
    solver_->SetLBSFissionSource(phi_temp_, false);
    Scale(q_moments_local_, 1.0 / lambda_k);

    CopyOnlyPhi0(q_moments_local_, Sf0_aux_);

    // Inner iterations seem extremely wasteful. Set this to 1 iteration here for further
    // investigation.
    for (int i = 0; i < 1; ++i)
    {
      SetLBSScatterSourcePhi0(epsilon_k_, false, true);

      CopyOnlyPhi0(q_moments_local_, Ss0_);

      // Solve the diffusion system
      diffusion_solver_->Assemble_b(Ss0_ + Sf0_aux_ + Ss0_res_ - Sf0_ell_);
      diffusion_solver_->Solve(epsilon_kp1_, true);

      epsilon_k_ = epsilon_kp1_;
    }

    ProjectBackPhi0(epsilon_kp1_ + phi0_star_, phi_old_local_);

    const double production_kp1 = ComputeFissionProduction(do_problem_, phi_old_local_);

    lambda_kp1 = production_kp1 / (production_k / lambda_k);

    const double lambda_change = std::fabs(lambda_kp1 / lambda_k - 1.0);
    if (verbose_ >= 1)
      log.Log() << "PISCDSA iteration " << k << " lambda " << lambda_kp1 << " lambda change "
                << lambda_change;

    if (lambda_change < pi_k_tol_)
      break;

    lambda_k = lambda_kp1;
    epsilon_k_ = epsilon_kp1_;
    production_k = production_kp1;
  } // acceleration loop

  ProjectBackPhi0(epsilon_kp1_ + phi0_star_, phi_new_local_);
  LBSVecOps::GSScopedCopyPrimarySTLvectors(
    do_problem_, front_gs_, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);

  return lambda_kp1;
}
} // namespace opensn
