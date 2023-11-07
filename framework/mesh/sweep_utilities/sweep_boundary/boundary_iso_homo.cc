#include "framework/mesh/sweep_utilities/sweep_boundary/boundary_iso_homo.h"

namespace chi_mesh::sweep_management
{

double*
BoundaryIsotropicHomogenous::HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                                      unsigned int face_num,
                                                      unsigned int fi,
                                                      unsigned int angle_num,
                                                      int group_num,
                                                      size_t gs_ss_begin)
{
  return &boundary_flux[group_num];
}

} // namespace chi_mesh::sweep_management
