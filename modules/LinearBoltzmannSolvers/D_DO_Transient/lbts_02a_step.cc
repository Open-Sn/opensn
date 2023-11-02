#include "lbts_transient_solver.h"

#include "LinearBoltzmannSolvers/B_DO_Solver/IterativeMethods/sweep_wgs_context.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Transient solver timestep routine.*/
void lbs::DiscOrdTransientSolver::Step()
{
  if (transient_options_.verbosity_level >= 2)
    chi::log.Log() << TextName() << " Stepping with dt " << dt_;

  phi_old_local_ = phi_prev_local_;

  for (auto& groupset : groupsets_)
  {
    //======================================== Converge the scattering source
    //                                         with a fixed fission source
    //                                         and temporal source
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    auto sweep_chunk = SetTransientSweepChunk(groupset);

    auto sweep_wgs_context_ptr =
    std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
      *this, groupset,
        active_set_source_function_,
        APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,  //lhs_scope
        APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
        APPLY_AGS_FISSION_SOURCES,                              //rhs_scope
        options_.verbose_inner_iterations,
        sweep_chunk);

    WGSLinearSolver<Mat,Vec,KSP> solver(sweep_wgs_context_ptr);
    solver.Setup();
    solver.Solve();

    Chi::mpi.Barrier();
  }

  //======================================== Compute t^{n+1} value
  {
    const auto& BackwardEuler = chi_math::SteppingMethod::IMPLICIT_EULER;
    const auto& CrankNicolson = chi_math::SteppingMethod::CRANK_NICOLSON;

    double theta;
    if      (method == BackwardEuler) theta = 1.0;
    else if (method == CrankNicolson) theta = 0.5;
    else                              theta = 0.7;
    const double inv_theta = 1.0/theta;

    auto& phi = phi_new_local_;
    const auto& phi_prev = phi_prev_local_;
    for (size_t i = 0; i < phi.size(); ++i)
      phi[i] = inv_theta*(phi[i] + (theta-1.0) * phi_prev[i]);

    if (options_.use_precursors)
      StepPrecursors();
  }

  const double FR_new = ComputeFissionProduction(phi_new_local_);

  //============================================= Print end of timestep
  if (transient_options_.verbosity_level >= 1)
  {
    char buff[200];
    snprintf(buff, 200, " dt=%.1e time=%10.4g FR=%12.6g", dt_, time_ + dt_, FR_new);
    chi::log.Log() << TextName() << buff;
  }

  UpdateFieldFunctions();
}