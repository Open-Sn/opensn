#include "sweep_boundaries.h"

// ###################################################################
/**Returns a pointer to a heterogeneous flux storage location.*/
double*
chi_mesh::sweep_management::BoundaryVaccuum::HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                                                      unsigned int face_num,
                                                                      unsigned int fi,
                                                                      unsigned int angle_num,
                                                                      int group_num,
                                                                      size_t gs_ss_begin)
{
  return &boundary_flux_[group_num];
}