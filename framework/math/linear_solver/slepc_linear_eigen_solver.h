// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_eigen_solver.h"
#include "framework/math/linear_solver/linear_solver_context.h"
#include <string>
#include <utility>
#include <memory>
#include <slepceps.h>

namespace opensn
{

class SLEPcLinearEigenSolver : public LinearEigenSolver
{
public:
  struct ToleranceOptions
  {
    double residual_absolute = 1.0e-6;
    int maximum_iterations = 100;
  } tolerance_options;

  SLEPcLinearEigenSolver(IterativeMethod method, std::shared_ptr<LinearEigenContext> context_ptr);

  virtual ~SLEPcLinearEigenSolver();

  ToleranceOptions& GetToleranceOptions() { return tolerance_options; }

  std::string GetEPSType() { return eps_type_; }

  void SetEPSType(std::string& type)
  {
    eps_type_ = type;
    if (type == "power")
      method_ = IterativeMethod::POWER;
    else if (type == "subspace")
      method_ = IterativeMethod::SUBSPACE;
    else if (type == "arnolidi")
      method_ = IterativeMethod::ARNOLDI;
    else if (type == "krylovschur")
      method_ = IterativeMethod::KRYLOV_SCHUR;
    else
      throw std::runtime_error("Invalid iterative method type for linear eigen solver");
  }

  void ApplyToleranceOptions();

  /// Set up the linaer solver
  virtual void Setup();

  /// Solve the system
  virtual void Solve();

protected:
  bool IsSystemSet() const { return system_set_; }
  virtual void PreSetupCallback();
  virtual void SetOptions();
  virtual void SetSolverContext();
  virtual void SetConvergenceTest();
  virtual void SetMonitor();
  virtual void SetPreconditioner();
  virtual void SetSystemSize() = 0;
  virtual void SetSystem() = 0;
  virtual void PostSetupCallback();
  virtual void PreSolveCallback();
  virtual void SetInitialGuess() = 0;
  virtual void PostSolveCallback();

  Mat A_;
  Vec x_;
  EPS eps_;
  std::string eps_type_;
  int64_t num_local_dofs_;
  int64_t num_global_dofs_;

private:
  bool system_set_;
};

} // namespace opensn
