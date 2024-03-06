#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/sweep_boundary.h"
#include <vector>
#include <limits>

namespace opensn
{
namespace lbs
{

/**
 * Reflective boundary condition.
 */
class ReflectingBoundary : public SweepBoundary
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
  std::vector<AngVec> boundary_flux_;
  std::vector<AngVec> boundary_flux_old_;

  std::vector<int> reflected_anglenum_;
  std::vector<std::vector<bool>> angle_readyflags_;

public:
  ReflectingBoundary(size_t num_groups,
                     const Normal& normal,
                     CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::REFLECTING, num_groups, coord_type), normal_(normal)
  {
  }

  const Vector3& Normal() const { return normal_; }

  bool IsOpposingReflected() const { return opposing_reflected_; }

  void SetOpposingReflected(bool value) { opposing_reflected_ = value; }

  std::vector<AngVec>& GetBoundaryFluxNew() { return boundary_flux_; }

  std::vector<AngVec>& GetBoundaryFluxOld() { return boundary_flux_old_; }

  std::vector<int>& GetReflectedAngleIndexMap() { return reflected_anglenum_; }

  std::vector<std::vector<bool>>& GetAngleReadyFlags() { return angle_readyflags_; }

  double* PsiIncoming(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int group_num,
                      size_t gs_ss_begin) override;

  double* PsiOutgoing(uint64_t cell_local_id,
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

} // namespace lbs
} // namespace opensn
