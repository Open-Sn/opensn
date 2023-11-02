#include "lbts_transient_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Constructor.*/
lbs::DiscOrdTransientSolver::DiscOrdTransientSolver(const std::string &in_text_name) :
  lbs::DiscOrdKEigenvalueSolver(in_text_name)
{
  chi::log.Log() << TextName() << " created.";
}

/**Destructor*/
lbs::DiscOrdTransientSolver::~DiscOrdTransientSolver()
{

}