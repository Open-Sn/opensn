#include "modules/LinearBoltzmannSolvers/D_DO_Transient/lbts_transient_solver.h"
#if 0
#include "modules/LinearBoltzmannSolvers/D_DO_Transient/SweepChunks/lbts_sweepchunk_pwl.h"

//###################################################################
/**Sets up the sweek chunk for the given discretization method.*/
std::shared_ptr<SweepChunk>
lbs::DiscOrdTransientSolver::SetTransientSweepChunk(LBSGroupset& groupset)
{
  double theta;
  if (method == chi_math::SteppingMethod::IMPLICIT_EULER) theta = 1.0;
  else
    theta = 0.5;

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk =
    std::make_shared<SweepChunkPWLTransientTheta>(grid_ptr_,             // Spatial grid of cells
                                                  *discretization_,      // Spatial discretization
                                                  unit_cell_matrices_,   // Unit cell matrices
                                                  cell_transport_views_, // Cell transport views
                                                  phi_new_local_,        // Destination phi
                                                  psi_new_local_[groupset.id_], // Destination psi

                                                  psi_prev_local_[groupset.id_],
                                                  theta,
                                                  dt_,

                                                  q_moments_local_, // Source moments
                                                  groupset,         // Reference groupset
                                                  matid_to_xs_map_, // Material cross-sections
                                                  num_moments_,
                                                  max_cell_dof_count_);

  return sweep_chunk;
}
#endif
