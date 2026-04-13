// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/cell/cell.h"
#include <vector>
#include <limits>
namespace opensn
{

/// Reflective boundary condition.
class ReflectingBoundary : public SweepBoundary
{
public:
  ReflectingBoundary(BoundaryBank& bank,
                     std::uint64_t bid,
                     const std::shared_ptr<MeshContinuum>& grid,
                     const std::vector<LBSGroupset>& groupsets,
                     const Vector3& normal,
                     CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN);

  const Vector3& GetNormal() const { return normal_; }

  const Vector3* GetNormalForReflection() const override { return &normal_; }

  bool HasDelayedAngularFlux() const override { return true; }

  /**
   * Get the list of anglesets that depend on a given angleset.
   * If the boundary is not opposing reflecting, this method extracts all anglesets in the angle
   * aggregation that can only begin sweeping after the given angle set has completed its sweep.
   * \param following_angle_sets Output set to which the dependent anglesets will be added.
   * \param angle_agg Angle aggregation containing all anglesets.
   * \param angleset Angleset for which dependent anglesets are sought.
   */
  void GetFollowingAngleSets(int groupset_id,
                             std::set<AngleSet*>& following_angle_sets,
                             const AngleAggregation& angle_agg,
                             const AngleSet& angleset) override;

  void SetOpposingReflected(
    std::uint64_t boundary_id,
    const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries) override;

  void ZeroOpposingDelayedAngularFluxOld() override;

  size_t CountDelayedAngularDOFsNew() const override;

  size_t CountDelayedAngularDOFsOld() const override;

  void AppendNewDelayedAngularDOFsToVector(std::vector<double>& output) const override;

  void AppendOldDelayedAngularDOFsToVector(std::vector<double>& output) const override;

  void AppendNewDelayedAngularDOFsToArray(int64_t& index, double* buffer) const override;

  void AppendOldDelayedAngularDOFsToArray(int64_t& index, double* buffer) const override;

  void SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* buffer) override;

  void SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* buffer) override;

  void SetNewDelayedAngularDOFsFromVector(const std::vector<double>& values,
                                          size_t& index) override;

  void SetOldDelayedAngularDOFsFromVector(const std::vector<double>& values,
                                          size_t& index) override;

  void CopyDelayedAngularFluxOldToNew() override;

  void CopyDelayedAngularFluxNewToOld() override;

  void GetReflectedMap(int groupset_id,
                       std::vector<std::vector<std::uint32_t>>& reflected_maps) override
  {
    reflected_maps.push_back(extra_data_[groupset_id].reflected_anglenum);
  }

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id,
                      unsigned int group_idx) override;

  double* PsiOutgoing(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id) override;

  std::uint64_t
  GetOffsetToAngleset(const FaceNode& face_node, AngleSet& anglset, bool is_outgoing) override;

  void InitializeReflectingMap(const std::vector<LBSGroupset>& groupsets) override;
  void InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets) override;

  /// Stride for computing psi pointer.
  struct Stride
  {
    std::uint64_t facenode_stride = 0;
    std::uint64_t angle_stride = 0;
  };

protected:
  const Vector3 normal_;
  CoordinateSystemType coord_type_;
  bool opposing_reflected_ = false;

  /// Map from face node to internal node-index.
  std::map<FaceNode, std::uint64_t> facenode_to_index_;

  /// Additional data specific to reflecting boundary for each groupset.
  struct ExtraData
  {
    /// Map from angle direction number in a quadrature to internal angle index
    std::vector<std::uint64_t> map_dirnum;
    /// Map from angle direction number in a quadrature to reflected angle index.
    std::vector<std::uint32_t> reflected_anglenum;
    /**
     * Number of angles recorded in the boundary, which is also the stride of the face node.
     * If the boundary is opposing reflected, all angles are recorded. Otherwise, only outgoing ones
     * are counted.
     */
    std::uint64_t node_stride = 0;
    /// Offset separating current and old boundary flux.
    std::uint64_t old_stride = 0;
  };
  /// List of data per groupset.
  std::vector<ExtraData> extra_data_;

private:
  template <bool ApplyOnOldFlux, typename Fn>
  void ForEachGroupset(Fn&& fn);

  template <bool ApplyOnOldFlux, typename Fn>
  void ForEachGroupset(Fn&& fn) const;

  template <bool ApplyOnOldFlux, typename Fn>
  void ForEachDelayedAngularFlux(Fn&& fn);

  template <bool ApplyOnOldFlux, typename Fn>
  void ForEachDelayedAngularFlux(Fn&& fn) const;

private:
  static constexpr double epsilon_ = 1.0e-8;
};

} // namespace opensn
