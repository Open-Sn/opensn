// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/lbts_transient_solver.h"
#if 0
#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/sweep_chunks/lbts_sweepchunk_pwl.h"

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
std::shared_ptr<SweepChunk>
lbs::DiscOrdTransientSolver::SetTransientSweepChunk(LBSGroupset& groupset)
{
  double theta;
  if (method == opensn::SteppingMethod::IMPLICIT_EULER) theta = 1.0;
  else
    theta = 0.5;

  // Setting up required sweep chunks
  auto sweep_chunk =
    std::make_shared<SweepChunkPWLTransientTheta>(grid_ptr_,
                                                  *discretization_,
                                                  unit_cell_matrices_,
                                                  cell_transport_views_,
                                                  phi_new_local_,
                                                  psi_new_local_[groupset.id_],

                                                  psi_prev_local_[groupset.id_],
                                                  theta,
                                                  dt_,

                                                  q_moments_local_,
                                                  groupset,
                                                  matid_to_xs_map_,
                                                  num_moments_,
                                                  max_cell_dof_count_);

  return sweep_chunk;
}
#endif
