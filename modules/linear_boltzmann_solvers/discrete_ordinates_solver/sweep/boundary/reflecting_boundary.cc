#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/reflecting_boundary.h"
#include "framework/logging/log.h"
#include "caliper/cali.h"

namespace opensn
{
namespace lbs
{

double*
ReflectingBoundary::PsiIncoming(uint64_t cell_local_id,
                                unsigned int face_num,
                                unsigned int fi,
                                unsigned int angle_num,
                                int group_num,
                                size_t gs_ss_begin)
{
  double* psi = nullptr;

  int reflected_angle_num = reflected_anglenum_[angle_num];

  if (opposing_reflected_)
  {
    psi = &boundary_flux_old_[reflected_angle_num][cell_local_id][face_num][fi][gs_ss_begin];
  }
  else
  {
    psi = &boundary_flux_[reflected_angle_num][cell_local_id][face_num][fi][gs_ss_begin];
  }

  return psi;
}

double*
ReflectingBoundary::PsiOutgoing(uint64_t cell_local_id,
                                unsigned int face_num,
                                unsigned int fi,
                                unsigned int angle_num,
                                size_t gs_ss_begin)
{
  return &boundary_flux_[angle_num][cell_local_id][face_num][fi][gs_ss_begin];
}

void
ReflectingBoundary::UpdateAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss)
{
  for (const size_t n : angles)
    angle_readyflags_[reflected_anglenum_[n]][gs_ss] = true;
}

bool
ReflectingBoundary::CheckAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss)
{
  if (opposing_reflected_)
    return true;
  bool ready_flag = true;
  for (auto& n : angles)
    if (not boundary_flux_[reflected_anglenum_[n]].empty())
      if (not angle_readyflags_[n][gs_ss])
        return false;

  return ready_flag;
}

void
ReflectingBoundary::ResetAnglesReadyStatus()
{
  boundary_flux_old_ = boundary_flux_;

  for (auto& flags : angle_readyflags_)
    for (int gs_ss = 0; gs_ss < flags.size(); ++gs_ss)
      flags[gs_ss] = false;
}

} // namespace lbs
} // namespace opensn
