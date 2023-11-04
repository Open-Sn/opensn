#include "boundary_reflecting.h"
#include "chi_log.h"

namespace chi_mesh::sweep_management
{

double*
BoundaryReflecting::HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                             unsigned int face_num,
                                             unsigned int fi,
                                             unsigned int angle_num,
                                             int group_num,
                                             size_t gs_ss_begin)
{
  double* Psi;

  int reflected_angle_num = reflected_anglenum_[angle_num];

  if (opposing_reflected_)
  {
    Psi = &hetero_boundary_flux_old_[reflected_angle_num][cell_local_id][face_num][fi][gs_ss_begin];
  }
  else
  {
    Psi = &hetero_boundary_flux_[reflected_angle_num][cell_local_id][face_num][fi][gs_ss_begin];
  }

  return Psi;
}

double*
BoundaryReflecting::HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                                             unsigned int face_num,
                                             unsigned int fi,
                                             unsigned int angle_num,
                                             size_t gs_ss_begin)
{
  return &hetero_boundary_flux_[angle_num][cell_local_id][face_num][fi][gs_ss_begin];
}

void
BoundaryReflecting::UpdateAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss)
{
  for (const size_t n : angles)
    angle_readyflags_[reflected_anglenum_[n]][gs_ss] = true;
}

bool
BoundaryReflecting::CheckAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss)
{
  if (opposing_reflected_) return true;
  bool ready_flag = true;
  for (auto& n : angles)
    if (!hetero_boundary_flux_[reflected_anglenum_[n]].empty())
      if (not angle_readyflags_[n][gs_ss]) return false;

  return ready_flag;
}

void
BoundaryReflecting::ResetAnglesReadyStatus()
{
  hetero_boundary_flux_old_ = hetero_boundary_flux_;

  for (auto& flags : angle_readyflags_)
    for (int gs_ss = 0; gs_ss < flags.size(); ++gs_ss)
      flags[gs_ss] = false;
}

} // namespace chi_mesh::sweep_management
