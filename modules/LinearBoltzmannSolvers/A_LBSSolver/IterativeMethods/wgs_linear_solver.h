#ifndef CHITECH_WGSLINSOLVEBASE_H
#define CHITECH_WGSLINSOLVEBASE_H

#include "math/LinearSolver/linear_solver.h"

#include "wgs_context.h"

#include <memory>
#include <vector>
#include <functional>

namespace lbs
{

//################################################################### Class def
/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
template<class MatType, class VecType, class SolverType>
class WGSLinearSolver : public
                        chi_math::LinearSolver<MatType,VecType,SolverType>
{
protected:
  std::vector<double> saved_q_moments_local_;

public:
  typedef std::shared_ptr<WGSContext<MatType,VecType,SolverType>> WGSContextPtr;

  /**Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.*/
  explicit WGSLinearSolver(WGSContextPtr gs_context_ptr) :
    chi_math::LinearSolver<MatType,VecType,SolverType>
      (IterativeMethodPETScName(gs_context_ptr->groupset_.iterative_method_),
       gs_context_ptr)
  {
    auto& groupset = gs_context_ptr->groupset_;
    auto& solver_tol_options = this->ToleranceOptions();
    solver_tol_options.residual_absolute  = groupset.residual_tolerance_;
    solver_tol_options.maximum_iterations = groupset.max_iterations_;
    solver_tol_options.gmres_restart_interval = groupset.gmres_restart_intvl_;
  }

protected:
  void PreSetupCallback() override;         //Customized via context
  /*virtual void SetOptions();*/
  /*virtual void SetSolverContext();*/      //Generic
  void SetConvergenceTest() override;       //Generic
  /*virtual void SetMonitor();*/

  virtual void SetSystemSize() override;    //Customized via context
  virtual void SetSystem() override;        //Generic

  void SetPreconditioner() override;        //Customized via context

  void PostSetupCallback() override;        //Customized via context
public:
  /*virtual void Setup();*/

protected:
  void PreSolveCallback() override;         //Customized via context
  void SetRHS() override;                   //Generic + with context elements
  void SetInitialGuess() override;          //Generic
  void PostSolveCallback() override;        //Generic + with context elements
public:
  /*virtual void Solve();*/

  virtual ~WGSLinearSolver() override;
};

}//namespace lbs

#endif //CHITECH_WGSLINSOLVEBASE_H
