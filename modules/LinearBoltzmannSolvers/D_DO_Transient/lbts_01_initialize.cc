#include "lbts_transient_solver.h"

#include "A_LBSSolver/SourceFunctions/transient_source_function.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Transient solver initialize routine.*/
void lbs::DiscOrdTransientSolver::Initialize()
{
  chi::log.Log() << "Initializing " << TextName() << ".";
  options_.save_angular_flux = true;
  DiscOrdKEigenvalueSolver::Initialize();
  DiscOrdKEigenvalueSolver::Execute();

  if (transient_options_.verbosity_level >= 1)
  {
    const double FR = ComputeFissionRate(phi_new_local_);
    char buff[200];
    snprintf(buff,200, " Initial Fission Rate FR=%12.6g", FR);
    chi::log.Log() << TextName() << buff;
  }

  //======================================== Compute auxiliary vectors
  fission_rate_local_.resize(grid_ptr_->local_cells.size(), 0.0);
  phi_prev_local_ = phi_old_local_;
  precursor_prev_local_ = precursor_new_local_;
  psi_prev_local_ = psi_new_local_;

  if (transient_options_.verbosity_level >= 0)
  {
    const double beta = ComputeBeta();
    char buff[200];
    snprintf(buff,200, " Beta=%.2f [pcm] reactivity=%.3f [$]",
            beta*1e5, (1.0- 1.0 / GetKeff()) / beta);
    chi::log.Log() << TextName() << buff;
  }

  //================================================== Initialize source func
  auto src_function =
    std::make_shared<TransientSourceFunction>(*this, this->dt_, this->method);

  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);
}