#include "lbs_curvilinear_solver.h"

#include "SweepChunks/lbs_curvilinear_sweepchunk_pwl.h"

namespace lbs
{

std::shared_ptr<chi_mesh::sweep_management::SweepChunk>
lbs::DiscreteOrdinatesCurvilinearSolver::SetSweepChunk(
  lbs::LBSGroupset& groupset)
{
  auto sweep_chunk =
    std::make_shared<SweepChunkPWLRZ>(*grid_ptr_,
                                      *discretization_,
                                      unit_cell_matrices_,
                                      secondary_unit_cell_matrices_,
                                      cell_transport_views_,
                                      phi_new_local_,
                                      psi_new_local_[groupset.id_],
                                      q_moments_local_,
                                      groupset,
                                      matid_to_xs_map_,
                                      num_moments_,
                                      max_cell_dof_count_);

  return sweep_chunk;
}

}