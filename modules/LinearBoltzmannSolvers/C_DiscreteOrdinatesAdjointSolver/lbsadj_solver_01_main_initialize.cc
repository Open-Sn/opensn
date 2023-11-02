#include "lbsadj_solver.h"

#include "A_LBSSolver/SourceFunctions/adjoint_src_function.h"

#include "chi_runtime.h"
#include "chi_log.h"

void
lbs::DiscreteOrdinatesAdjointSolver::Initialize()
{
  LBSSolver::Initialize();

  MakeAdjointXSs();
  InitQOIs();

  //================================================== Initialize source func
  auto src_function = std::make_shared<AdjointSourceFunction>(*this);

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

  InitializeSolverSchemes(); // j
  source_event_tag_ = Chi::log.GetRepeatingEventTag("Set Source");
}