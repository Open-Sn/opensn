#include "pi_keigen_scdsa.h"

#include "A_LBSSolver/Acceleration/diffusion_mip.h"
#include "A_LBSSolver/Acceleration/diffusion_PWLC.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "math/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs
{

void
XXPowerIterationKEigenSCDSA::Initialize()
{
  XXPowerIterationKEigen::Initialize();

  //=========================================== Make UnknownManager
  const size_t num_gs_groups = front_gs_.groups_.size();
  chi_math::UnknownManager uk_man;
  uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, num_gs_groups);

  //=========================================== Make boundary conditions
  auto bcs = acceleration::TranslateBCs(lbs_solver_.SweepBoundaries(),
                                        /*vaccum_bcs_are_dirichlet=*/true);

  //=========================================== Make xs map
  auto matid_2_mgxs_map = acceleration::PackGroupsetXS(
    lbs_solver_.GetMatID2XSMap(), front_gs_.groups_.front().id_, front_gs_.groups_.back().id_);

  //=========================================== Create solver
  const auto& sdm = lbs_solver_.SpatialDiscretization();
  const auto& unit_cell_matrices = lbs_solver_.GetUnitCellMatrices();

  if (diffusion_solver_sdm_ == "pwld")
    diffusion_solver_ =
      std::make_shared<acceleration::DiffusionMIPSolver>(std::string(TextName() + "_WGDSA"),
                                                         sdm,
                                                         uk_man,
                                                         bcs,
                                                         matid_2_mgxs_map,
                                                         unit_cell_matrices,
                                                         true); // verbosity
  else
  {
    continuous_sdm_ptr_ =
      chi_math::spatial_discretization::PieceWiseLinearContinuous::New(sdm.Grid());
    diffusion_solver_ =
      std::make_shared<acceleration::DiffusionPWLCSolver>(std::string(TextName() + "_WGDSA"),
                                                          *continuous_sdm_ptr_,
                                                          uk_man,
                                                          bcs,
                                                          matid_2_mgxs_map,
                                                          unit_cell_matrices,
                                                          true); // verbosity
    requires_ghosts_ = true;
    lbs_pwld_ghost_info_ =
      MakePWLDVecGhostCommInfo(lbs_solver_.SpatialDiscretization(), lbs_solver_.UnknownManager());

    const auto& cfem_sdm = *continuous_sdm_ptr_;
    const auto ghost_dof_ids = cfem_sdm.GetGhostDOFIndices(lbs_solver_.UnknownManager());
  }

  {
    auto& ds = diffusion_solver_;

    ds->options.residual_tolerance = diff_accel_diffusion_l_abs_tol_;
    ds->options.max_iters = diff_accel_diffusion_max_iters_;
    ds->options.verbose = diff_accel_diffusion_verbose_;
    ds->options.additional_options_string = diff_accel_diffusion_petsc_options_;
  }

  Chi::log.Log() << "Initializing diffusion solver";
  diffusion_solver_->Initialize();
  Chi::mpi.Barrier();
  Chi::log.Log() << "Done Initializing diffusion solver";

  Chi::log.Log() << "Assembling A and b";
  std::vector<double> dummy_rhs;
  if (diffusion_solver_sdm_ == "pwld") dummy_rhs.assign(sdm.GetNumLocalDOFs(uk_man), 0.0);
  else
    dummy_rhs.assign(continuous_sdm_ptr_->GetNumLocalAndGhostDOFs(uk_man), 0.0);

  diffusion_solver_->AssembleAand_b(dummy_rhs);
  Chi::log.Log() << "Done Assembling A and b";
}

} // namespace lbs