// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/source_function.h"

namespace opensn
{

/**
 * A transient source function needs to adjust the DelayedFission routine to properly fit with
 * the current timestepping method and timestep.
 */
class TransientSourceFunction : public SourceFunction
{
public:
  /**
   * Constructor for the transient source function. The only difference as compared to a steady
   * source function is the treatment of delayed fission.
   */
  explicit TransientSourceFunction(const LBSProblem& lbs_problem);

  double DelayedFission(const PrecursorList& precursors,
                        const double& rho,
                        const std::vector<double>& nu_delayed_sigma_f,
                        const double* phi) const override;

private:
  const LBSProblem& lbs_problem_;
};

} // namespace opensn
