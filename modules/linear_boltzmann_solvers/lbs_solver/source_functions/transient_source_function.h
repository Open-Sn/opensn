// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/source_functions/source_function.h"
#include "framework/math/math_time_stepping.h"

namespace opensn
{

/**
 * A transient source function needs to adjust the AddDelayedFission routine to properly fit with
 * the current timestepping method and timestep.
 */
class TransientSourceFunction : public SourceFunction
{
private:
  double& dt_;
  SteppingMethod& method_;

public:
  /**
   * Constructor for the transient source function. The only difference as compared to a steady
   * source function is the treatment of delayed fission.
   */
  TransientSourceFunction(const LBSSolver& lbs_solver, double& ref_dt, SteppingMethod& method);

  double AddDelayedFission(const PrecursorList& precursors,
                           const double& rho,
                           const std::vector<double>& nu_delayed_sigma_f,
                           const double* phi) const override;
};

} // namespace opensn
