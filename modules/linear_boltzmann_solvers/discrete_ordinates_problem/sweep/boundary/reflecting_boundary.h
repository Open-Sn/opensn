// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include <vector>
#include <limits>

namespace opensn
{

/// Reflective boundary condition.
class ReflectingBoundary : public SweepBoundary
{
protected:
  const Vector3 normal_;
  bool opposing_reflected_ = false;

  // Populated by angle aggregation
  // AngularFluxData indices are: angle, cell per angle, faces per cell, DOFs per face, groups per
  // DOF
  using AngularFluxData = std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>;
  AngularFluxData boundary_flux_;
  AngularFluxData boundary_flux_old_;

  std::vector<int> reflected_anglenum_;
  std::vector<bool> angle_readyflags_;

public:
  ReflectingBoundary(size_t num_groups,
                     const Vector3& normal,
                     CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(LBSBoundaryType::REFLECTING, num_groups, coord_type), normal_(normal)
  {
  }

  const Vector3& GetNormal() const { return normal_; }

  bool IsOpposingReflected() const { return opposing_reflected_; }

  void SetOpposingReflected(bool value) { opposing_reflected_ = value; }

  AngularFluxData& GetBoundaryFluxNew() { return boundary_flux_; }

  AngularFluxData& GetBoundaryFluxOld() { return boundary_flux_old_; }

  std::vector<int>& GetReflectedAngleIndexMap() { return reflected_anglenum_; }

  std::vector<bool>& GetAngleReadyFlags() { return angle_readyflags_; }

  double* PsiIncoming(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int group_num) override;

  double* PsiOutgoing(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num) override;

  void UpdateAnglesReadyStatus(const std::vector<size_t>& angles) override;

  bool CheckAnglesReadyStatus(const std::vector<size_t>& angles) override;

  /// Resets angle ready flags to false.
  void ResetAnglesReadyStatus();
};

} // namespace opensn
