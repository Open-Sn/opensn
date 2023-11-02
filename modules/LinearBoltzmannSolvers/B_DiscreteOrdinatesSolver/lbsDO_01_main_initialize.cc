#include "lbs_discrete_ordinates_solver.h"

#include "B_DiscreteOrdinatesSolver/IterativeMethods/sweep_wgs_context.h"
#include "A_LBSSolver/IterativeMethods/wgs_linear_solver.h"
#include "A_LBSSolver/SourceFunctions/source_function.h"

#include "chi_runtime.h"
#include "chi_log.h"

typedef chi_mesh::sweep_management::SweepChunk SweepChunk;

//###################################################################
/** Initialize the solver.*/
void lbs::DiscreteOrdinatesSolver::Initialize()
{
  LBSSolver::Initialize();

  auto src_function = std::make_shared<SourceFunction>(*this);

  // Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);

  //================================================== Initialize groupsets for
  //                                                   sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  InitializeSolverSchemes();           //j
  source_event_tag_ = Chi::log.GetRepeatingEventTag("Set Source");
}

/**Initializes Within-GroupSet solvers.*/
void lbs::DiscreteOrdinatesSolver::InitializeWGSSolvers()
{
  wgs_solvers_.clear(); //this is required
  for (auto& groupset : groupsets_)
  {
    std::shared_ptr<SweepChunk> sweep_chunk = SetSweepChunk(groupset);

    auto sweep_wgs_context_ptr =
    std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
      *this, groupset,
        active_set_source_function_,
        APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,  //lhs_scope
        APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
        APPLY_AGS_FISSION_SOURCES,                              //rhs_scope
        options_.verbose_inner_iterations,
        sweep_chunk);

    auto wgs_solver =
      std::make_shared<WGSLinearSolver<Mat,Vec,KSP>>(sweep_wgs_context_ptr);

    wgs_solvers_.push_back(wgs_solver);
  }//for groupset

}