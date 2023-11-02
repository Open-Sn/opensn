#include "pi_keigen_scdsa.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs
{

RegisterChiObject(lbs, XXPowerIterationKEigenSCDSA);

chi::InputParameters XXPowerIterationKEigenSCDSA::GetInputParameters()
{
  chi::InputParameters params =
    XXPowerIterationKEigen::GetInputParameters();

  params.SetGeneralDescription(
    "Generalized implementation of a k-Eigenvalue solver using Power "
    "Iteration and with SCDSA acceleration.");
  params.SetDocGroup("LBSExecutors");

  params.AddOptionalParameter(
    "accel_pi_max_its",
    50,
    "Maximum allowable iterations for the acceleration scheme's inner "
    "power iterations");

  params.AddOptionalParameter(
    "accel_pi_k_tol",
    1.0e-10,
    "K-eigenvalue tolerance for the acceleration scheme's inner "
    "power iterations");

  params.AddOptionalParameter(
    "accel_pi_verbose",
    false,
    "Flag, if set will result in verbose output from the acceleration "
    "scheme");

  params.AddOptionalParameter(
    "diff_accel_diffusion_l_abs_tol",
    1.0e-10,
    "Absolute residual tolerance to use for the diffusion accelerator");
  params.AddOptionalParameter(
    "diff_accel_diffusion_max_iters",
    100,
    "Maximum allowable iterations for the diffusion accelerator");
  params.AddOptionalParameter(
    "diff_accel_diffusion_verbose",
    false,
    "Flag, if set will enable verbose output of the diffusion accelerator");
  params.AddOptionalParameter(
    "diff_accel_diffusion_petsc_options",
    std::string("ssss"),
    "Additional PETSc options for the diffusion accelerator");

  params.AddOptionalParameter(
    "diff_accel_sdm",
    "pwld",
    "Spatial discretization to use for the diffusion solver");

  using namespace chi_data_types;
  params.ConstrainParameterRange("diff_accel_sdm",
                                 AllowableRangeList::New({"pwld", "pwlc"}));

  return params;
}

XXPowerIterationKEigenSCDSA::XXPowerIterationKEigenSCDSA(
  const chi::InputParameters& params)
  : XXPowerIterationKEigen(params),
    accel_pi_max_its_(params.GetParamValue<int>("accel_pi_max_its")),
    accel_pi_k_tol_(params.GetParamValue<double>("accel_pi_k_tol")),
    accel_pi_verbose_(params.GetParamValue<bool>("accel_pi_verbose")),
    diffusion_solver_sdm_(params.GetParamValue<std::string>("diff_accel_sdm")),
    diff_accel_diffusion_l_abs_tol_(
      params.GetParamValue<double>("diff_accel_diffusion_l_abs_tol")),
    diff_accel_diffusion_max_iters_(
      params.GetParamValue<int>("diff_accel_diffusion_max_iters")),
    diff_accel_diffusion_verbose_(
      params.GetParamValue<bool>("diff_accel_diffusion_verbose")),
    diff_accel_diffusion_petsc_options_(
      params.GetParamValue<std::string>("diff_accel_diffusion_petsc_options"))
{
  ////=========================================== Make UnknownManager
  // const size_t num_gs_groups = front_gs_.groups_.size();
  // chi_math::UnknownManager uk_man;
  // uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_gs_groups);
  //
  ////=========================================== Make boundary conditions
  // auto bcs = acceleration::TranslateBCs(lbs_solver_.SweepBoundaries(),
  //                                       /*vaccum_bcs_are_dirichlet=*/true);
  //
  ////=========================================== Make xs map
  // auto matid_2_mgxs_map =
  //   acceleration::PackGroupsetXS(lbs_solver_.GetMatID2XSMap(),
  //                                front_gs_.groups_.front().id_,
  //                                front_gs_.groups_.back().id_);
  //
  ////=========================================== Create solver
  // const auto& sdm = lbs_solver_.SpatialDiscretization();
  // const auto& unit_cell_matrices = lbs_solver_.GetUnitCellMatrices();
  //
  // if (diffusion_solver_sdm_ == "pwld")
  //   diffusion_solver_ = std::make_shared<acceleration::DiffusionMIPSolver>(
  //     std::string(TextName() + "_WGDSA"),
  //     sdm,
  //     uk_man,
  //     bcs,
  //     matid_2_mgxs_map,
  //     unit_cell_matrices,
  //     true); // verbosity
  // else
  //{
  //   continuous_sdm_ptr_ =
  //     chi_math::SpatialDiscretization_PWLC::New(sdm.ref_grid_);
  //   diffusion_solver_ = std::make_shared<acceleration::DiffusionPWLCSolver>(
  //     std::string(TextName() + "_WGDSA"),
  //     *continuous_sdm_ptr_,
  //     uk_man,
  //     bcs,
  //     matid_2_mgxs_map,
  //     unit_cell_matrices,
  //     true); // verbosity
  //   requires_ghosts_ = true;
  //   lbs_pwld_ghost_info_ = MakePWLDVecGhostCommInfo(
  //     lbs_solver_.SpatialDiscretization(), lbs_solver_.UnknownManager());
  //
  //   const auto& cfem_sdm = *continuous_sdm_ptr_;
  //   const auto ghost_dof_ids =
  //     cfem_sdm.GetGhostDOFIndices(lbs_solver_.UnknownManager());
  // }
  //
  //{
  //   typedef const std::string cstr;
  //   cstr l_abs_tol = "diff_accel_diffusion_l_abs_tol";
  //   cstr max_iters = "diff_accel_diffusion_max_iters";
  //   cstr verbose = "diff_accel_diffusion_verbose";
  //   cstr petsc_options = "diff_accel_diffusion_petsc_options";
  //
  //   auto& ds = diffusion_solver_;
  //
  //   ds->options.residual_tolerance = params.GetParamValue<double>(l_abs_tol);
  //   ds->options.max_iters = params.GetParamValue<int>(max_iters);
  //   ds->options.verbose = params.GetParamValue<bool>(verbose);
  //   if (not params.GetParamValue<std::string>(petsc_options).empty())
  //     ds->options.additional_options_string =
  //       params.GetParamValue<std::string>(petsc_options);
  // }
  //
  // chi::log.Log() << "Initializing diffusion solver";
  // diffusion_solver_->Initialize();
  // Chi::mpi.Barrier();
  // chi::log.Log() << "Done Initializing diffusion solver";
  //
  // chi::log.Log() << "Assembling A and b";
  // std::vector<double> dummy_rhs;
  // if (diffusion_solver_sdm_ == "pwld")
  //   dummy_rhs.assign(sdm.GetNumLocalDOFs(uk_man), 0.0);
  // else
  //   dummy_rhs.assign(continuous_sdm_ptr_->GetNumLocalAndGhostDOFs(uk_man),
  //   0.0);
  //
  // diffusion_solver_->AssembleAand_b(dummy_rhs);
  // chi::log.Log() << "Done Assembling A and b";
}

} // namespace lbs