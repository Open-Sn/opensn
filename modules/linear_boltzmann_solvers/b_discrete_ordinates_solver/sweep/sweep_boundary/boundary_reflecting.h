#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/sweep_boundary/sweep_boundary.h"
#include <vector>
#include <limits>

namespace opensn
{

/**
 * Reflective boundary condition.
 */
class BoundaryReflecting : public SweepBoundary
{
protected:
  const opensn::Normal normal_;
  bool opposing_reflected_ = false;

  /// Groups per DOF
  typedef std::vector<double> DOFVec;
  /// DOFs per face
  typedef std::vector<DOFVec> FaceVec;
  /// Faces per cell
  typedef std::vector<FaceVec> CellVec;
  /// Cell per angle
  typedef std::vector<CellVec> AngVec;

  // angle,cell,face,dof,group
  // Populated by angle aggregation
  std::vector<AngVec> hetero_boundary_flux_;
  std::vector<AngVec> hetero_boundary_flux_old_;

  std::vector<int> reflected_anglenum_;
  std::vector<std::vector<bool>> angle_readyflags_;

public:
  BoundaryReflecting(size_t in_num_groups,
                     const Normal& in_normal,
                     CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::REFLECTING, in_num_groups, coord_type), normal_(in_normal)
  {
  }

  const Vector3& Normal() const { return normal_; }
  bool IsOpposingReflected() const { return opposing_reflected_; }
  void SetOpposingReflected(bool value) { opposing_reflected_ = value; }

  std::vector<AngVec>& GetHeteroBoundaryFluxNew() { return hetero_boundary_flux_; }
  std::vector<AngVec>& GetHeteroBoundaryFluxOld() { return hetero_boundary_flux_old_; }

  std::vector<int>& GetReflectedAngleIndexMap() { return reflected_anglenum_; }
  std::vector<std::vector<bool>>& GetAngleReadyFlags() { return angle_readyflags_; }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;
  double* HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   size_t gs_ss_begin) override;

  void UpdateAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss) override;
  bool CheckAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss) override;

  /**
   * Resets angle ready flags to false.
   */
  void ResetAnglesReadyStatus();
};

} // namespace opensn
