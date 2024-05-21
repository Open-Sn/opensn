// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_context.h"

namespace opensn
{
namespace lbs
{

/**Linear Solver specialization for Within GroupSet (WGS) solves.*/
class AGSLinearSolver : public LinearSolver
{
public:
  /**Constructor.
   * \param iterative_method string Across Groupset iterative method.
   * \param ags_context_ptr Pointer Pointer to the context to use.
   * \param groupspan_first_id int First group index.
   * \param groupspan_last_id int Last group index.
   * \param verbose bool Flag to enable verbose output.*/
  AGSLinearSolver(std::string iterative_method,
                  std::shared_ptr<AGSContext> ags_context_ptr,
                  int groupspan_first_id,
                  int groupspan_last_id,
                  bool verbose = true)
    : LinearSolver(std::move(iterative_method), ags_context_ptr),
      groupspan_first_id_(groupspan_first_id),
      groupspan_last_id_(groupspan_last_id),
      verbose_(verbose)
  {
  }

  ~AGSLinearSolver() override;

  int GroupSpanFirstID() const { return groupspan_first_id_; }

  int GroupSpanLastID() const { return groupspan_last_id_; }

  bool IsVerbose() const { return verbose_; }

  void SetVerbosity(bool verbose_y_n) { verbose_ = verbose_y_n; }

  void Solve() override;

protected:
  void SetSystemSize() override;

  void SetSystem() override;

  void SetPreconditioner() override;

  void SetRHS() override;

  void SetInitialGuess() override;

  int groupspan_first_id_;

  int groupspan_last_id_;

  bool verbose_;
};

} // namespace lbs
} // namespace opensn
