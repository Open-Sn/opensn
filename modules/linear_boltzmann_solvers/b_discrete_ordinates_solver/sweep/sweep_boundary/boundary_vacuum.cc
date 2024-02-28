#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/sweep_boundary/boundary_vacuum.h"

namespace opensn
{
namespace lbs
{

double*
BoundaryVaccuum::HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                          unsigned int face_num,
                                          unsigned int fi,
                                          unsigned int angle_num,
                                          int group_num,
                                          size_t gs_ss_begin)
{
  return &boundary_flux_[group_num];
}

} // namespace lbs
} // namespace opensn
