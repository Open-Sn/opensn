// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_bank.h"
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace opensn
{

class AngleSet;
class AngleAggregation;
class FaceNode;
class MeshContinuum;
class LBSGroupset;

/// Base class for sweep related boundaries.
class SweepBoundary
{
public:
  explicit SweepBoundary(BoundaryBank& bank, LBSBoundaryType bndry_type)
    : bank_(bank), type_(bndry_type)
  {
    offset_.reserve(bank_.GetSize());
    offset_.resize(bank_.GetSize());
  }

  virtual ~SweepBoundary() = default;

  LBSBoundaryType GetType() const { return type_; }

  bool IsReflecting() const { return type_ == LBSBoundaryType::REFLECTING; }
  bool IsAngleDependent() const
  {
    return type_ == LBSBoundaryType::REFLECTING || type_ == LBSBoundaryType::ARBITRARY;
  }

  double GetEvaluationTime() const { return evaluation_time_; }

  void SetEvaluationTime(double time) { evaluation_time_ = time; }

  virtual bool HasDelayedAngularFlux() const { return false; }

  /**
   * Get the list of anglesets that depend on a given angleset.
   * This function can only be applied to reflecting boundary and AAHD anglesets.
   */
  virtual void GetFollowingAngleSets(int groupset_id,
                                     std::set<AngleSet*>& following_angle_sets,
                                     const AngleAggregation& angle_agg,
                                     const AngleSet& angleset)
  {
  }

  virtual void
  SetOpposingReflected(std::uint64_t boundary_id,
                       const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
  {
  }

  virtual void ZeroOpposingDelayedAngularFluxOld() {}

  virtual size_t CountDelayedAngularDOFsNew() const { return 0; }

  virtual size_t CountDelayedAngularDOFsOld() const { return 0; }

  virtual void AppendNewDelayedAngularDOFsToVector(std::vector<double>& output) const {}

  virtual void AppendOldDelayedAngularDOFsToVector(std::vector<double>& output) const {}

  virtual void AppendNewDelayedAngularDOFsToArray(int64_t& index, double* buffer) const {}

  virtual void AppendOldDelayedAngularDOFsToArray(int64_t& index, double* buffer) const {}

  virtual void SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* buffer) {}

  virtual void SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* buffer) {}

  virtual void SetNewDelayedAngularDOFsFromVector(const std::vector<double>& values, size_t& index)
  {
  }

  virtual void SetOldDelayedAngularDOFsFromVector(const std::vector<double>& values, size_t& index)
  {
  }

  virtual void CopyDelayedAngularFluxOldToNew() {}

  virtual void CopyDelayedAngularFluxNewToOld() {}

  virtual const Vector3* GetNormalForReflection() const { return nullptr; }

  virtual void GetReflectedMap(int groupset_id,
                               std::vector<std::vector<std::uint32_t>>& reflected_maps)
  {
  }

  /// Return a pointer to the location of the incoming flux.
  virtual double* PsiIncoming(std::uint32_t cell_local_id,
                              unsigned int face_num,
                              unsigned int fi,
                              unsigned int angle_num,
                              int groupset_id,
                              unsigned int group_idx);

  /// Return a pointer to the location of the outgoing flux.
  virtual double* PsiOutgoing(uint64_t cell_local_id,
                              unsigned int face_num,
                              unsigned int fi,
                              unsigned int angle_num,
                              int groupset_id);

  virtual void Setup(const std::shared_ptr<MeshContinuum>& grid,
                     const AngularQuadrature& quadrature)
  {
  }

  virtual void InitializeReflectingMap(const std::vector<LBSGroupset>& groupsets) {}

  /**
   * Perform additional setup required by angle-dependent boundaries.
   * Angle-dependent boundaries require an extra initialization step after flux data structures are
   * initialized. This phase cannot be executed in the constructor, which happens before the
   * initialization of flux data structures.
   */
  virtual void InitializeAngleDependent(const std::vector<LBSGroupset>& groupsets) {}

  /// Get pointer to the boundary flux section in the contigous bank for a given groupset.
  double* GetBoundaryFlux(int groupset_id, std::uint64_t extra_offset = 0)
  {
    auto& data = bank_[groupset_id];
    return data.boundary_flux.data() + (offset_[groupset_id] + extra_offset) * data.groupset_size;
  }

  /// Get pointer to the boundary flux section in the contigous bank for a given groupset.
  const double* GetBoundaryFlux(int groupset_id, std::uint64_t extra_offset = 0) const
  {
    auto& data = bank_[groupset_id];
    return data.boundary_flux.data() + (offset_[groupset_id] + extra_offset) * data.groupset_size;
  }

  /// Get offset to the data region of a given angleset.
  virtual std::uint64_t
  GetOffsetToAngleset(const FaceNode& face_node, AngleSet& anglset, bool is_outgoing)
  {
    return std::numeric_limits<std::uint64_t>::max();
  }

  double* ZeroFlux(int groupset_id, unsigned int group_num)
  {
    auto& data = bank_[groupset_id];
    return data.boundary_flux.data() + group_num;
  }

protected:
  BoundaryBank& bank_;

  /**
   * Offset used to access boundary flux data for each groupset.
   * Boundary fluxes for a given boundary start at offset * groupset.GetNumGroups().
   */
  std::vector<std::uint64_t> offset_;

private:
  /// Boundary type.
  LBSBoundaryType type_;
  /// Time value passed to boundary functions
  double evaluation_time_ = 0.0;
};

/**
 * This boundary function class can be derived from to
 * provide a much more custom experience. This function
 * is called during Setup.
 */
class BoundaryFunction
{
public:
  /// Customized boundary function by calling an input routine.
  virtual std::vector<double>
  Evaluate(size_t cell_global_id,
           int cell_block_id,
           unsigned int face_index,
           unsigned int face_node_index,
           const Vector3& face_node_location,
           const Vector3& face_node_normal,
           const std::vector<int>& quadrature_angle_indices,
           const std::vector<Vector3>& quadrature_angle_vectors,
           const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
           const std::vector<int>& group_indices,
           double time) = 0;

  virtual ~BoundaryFunction() = default;
};

} // namespace opensn
