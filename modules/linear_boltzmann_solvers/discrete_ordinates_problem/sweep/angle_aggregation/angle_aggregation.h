// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/mesh/mesh/mesh.h"
#include <memory>
#include <map>
#include <set>
#include <vector>

namespace opensn
{

using DirIDs = std::vector<size_t>; ///< Direction-IDs
using UniqueSOGroupings = std::vector<DirIDs>;
using DirIDToSOMap = std::map<size_t, size_t>;

enum class AngleAggregationType
{
  UNDEFINED = 0,
  SINGLE = 1,
  POLAR = 2,
  AZIMUTHAL = 3,
  QUADRANT = 4,
  OCTANT = 5,
};

/**
 * Angle aggregation partitions the directions of a quadrature into angle sets. All directions in
 * an angle set share one SPDS (sweep ordering), so the set is swept together. The
 * AngleAggregationType determines the partition:
 *
 * - SINGLE: one direction per angle set.
 * - POLAR: directions with the same azimuthal angle in the same polar hemisphere.
 * - AZIMUTHAL: the inward or outward directions of a curvilinear direction set, in angular order.
 * - QUADRANT: directions in the same quadrant (same signs of Omega_x and Omega_y), on a 2D
 *   orthogonal mesh.
 * - OCTANT: directions in the same octant (same signs of Omega_x, Omega_y, and Omega_z), on a 3D
 *   orthogonal mesh.
 *
 * On an orthogonal mesh, directions in the same quadrant/octant have the same orientation
 * (incoming or outgoing) on every face.
 *
 * The angle sets of a groupset are stored in a flat list (GetAngleSetGroups()) that the sweep
 * scheduler iterates over.
 */
class AngleAggregation
{
public:
  AngleAggregation(const LBSGroupset& groupset,
                   std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                   std::shared_ptr<AngularQuadrature>& quadrature,
                   std::shared_ptr<Mesh>& grid,
                   bool csda_enabled);

  using iterator = std::vector<std::shared_ptr<AngleSet>>::iterator;
  using const_iterator = std::vector<std::shared_ptr<AngleSet>>::const_iterator;
  iterator begin() { return angle_set_groups_.begin(); }
  iterator end() { return angle_set_groups_.end(); }
  const_iterator begin() const { return angle_set_groups_.begin(); }
  const_iterator end() const { return angle_set_groups_.end(); }
  std::shared_ptr<AngleSet> operator[](std::size_t idx) const { return angle_set_groups_[idx]; }

  std::vector<std::shared_ptr<AngleSet>>& GetAngleSetGroups() { return angle_set_groups_; }
  const std::vector<std::shared_ptr<AngleSet>>& GetAngleSetGroups() const
  {
    return angle_set_groups_;
  }
  std::size_t GetNumAngleSets() const { return angle_set_groups_.size(); }

  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetSimBoundaries() const
  {
    return boundaries_;
  }
  std::shared_ptr<AngularQuadrature> GetQuadrature() const { return quadrature_; }
  CoordinateSystemType GetCoordinateSystem() const
  {
    return grid_ ? grid_->GetCoordinateSystem() : CoordinateSystemType::UNDEFINED;
  }

  int GetGroupsetID() const { return groupset_id_; }

  /// Resets all the outgoing intra-location and inter-location cyclic interfaces.
  void ZeroOutgoingDelayedPsi();

  /// Resets all the incoming intra-location and inter-location cyclic interfaces.
  void ZeroIncomingDelayedPsi();

  /// Resets dependencies between anglesets.
  void ResetAngleSetDependencies();

  /// Establishes dependencies between anglesets for reflecting boundaries and cylindrical geometry.
  void SetupAngleSetDependencies();

  /// Returns a map of following angle sets per angle set.
  const std::map<AngleSet*, std::set<AngleSet*>>& GetFollowingAngleSetsMap() const
  {
    return following_angle_sets_map_;
  }

  /// Builds the map from direction number to angleset pointer
  void BuildDirnumToAnglesetMap();

  /// Gets the angle set for a given angle index.
  AngleSet* GetAngleSetForAngleIndex(std::uint32_t angle_idx) const
  {
    return dir_to_angleset_.at(angle_idx);
  }

  /**
   * Returns a pair of numbers containing the number of delayed angular unknowns both locally and
   * globally, respectively.
   */
  std::pair<size_t, size_t> GetNumDelayedAngularDOFs();

  /// Assembles angular unknowns into the reference vector.
  void AppendNewDelayedAngularDOFsToArray(int64_t& index, double* x_ref);

  /// Assembles angular unknowns into the reference vector.
  void AppendOldDelayedAngularDOFsToArray(int64_t& index, double* x_ref);

  /// Assembles angular unknowns into the reference vector.
  void SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref);

  /// Assembles angular unknowns into the reference vector.
  void SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref);

  /// Gets the current values of the angular unknowns as an STL vector.
  std::vector<double> GetNewDelayedAngularDOFsAsSTLVector();

  /// Gets the current values of the angular unknowns as an STL vector.
  void SetNewDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector);

  /// Gets the current values of the angular unknowns as an STL vector.
  std::vector<double> GetOldDelayedAngularDOFsAsSTLVector();

  /// Gets the current values of the angular unknowns as an STL vector.
  void SetOldDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector);

  /// Copies the old delayed angular fluxes to the new.
  void SetDelayedPsiOld2New();

  /// Copies the new delayed angular fluxes to the old.
  void SetDelayedPsiNew2Old();

  /**
   * Adds an isotropic nodal correction to every new delayed angular flux of this groupset: the
   * lagged incoming fluxes of opposing reflecting boundaries and of broken cyclic dependencies
   * within and between locations.
   */
  void CorrectNewDelayedAngularFlux(const FLUDS::DelayedFluxCorrection& correction);

private:
  int groupset_id_;
  bool num_ang_unknowns_avail_;
  std::pair<size_t, size_t> number_angular_unknowns_;
  bool csda_enabled_;
  std::shared_ptr<Mesh> grid_;
  std::shared_ptr<AngularQuadrature> quadrature_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries_;
  std::vector<std::shared_ptr<AngleSet>> angle_set_groups_;
  std::map<std::uint32_t, AngleSet*> dir_to_angleset_;
  std::map<AngleSet*, std::set<AngleSet*>> following_angle_sets_map_;
};

} // namespace opensn
