// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set_group.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <memory>
#include <map>
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
};

/**
 * Angle aggregation has to cater for running the 8 corners of a 3D partitioning, the 4 corners of a
 * 2D partitioning (the latter 2 both being polar angle aggregation) as well as single angle
 * aggregation.
 *
 * At the most fundamental level this manifests as a number of angle indices that share a SPDS,
 * however SPDS do not have to be unique which allows the notion of polar angle sets. For single
 * angle aggregation each SPDS is associated with a single angle index. For polar angle aggregation
 * a single SPDS can have multiple indices associated with the same azimuthal angle but different
 * polar angles. We call this manifestation a "AngleSet".
 *
 * The octant based separation is achieved via the notion of "AngleSetGroup" which will group angle
 * sets for each quadrant or octant (depending on 2D or 3D).
 */
class AngleAggregation
{
private:
  size_t num_groups_;
  bool num_ang_unknowns_avail_;
  std::pair<size_t, size_t> number_angular_unknowns_;
  std::shared_ptr<MeshContinuum> grid_;
  std::shared_ptr<AngularQuadrature> quadrature_;
  std::map<uint64_t, std::shared_ptr<SweepBoundary>> boundaries_;

public:
  AngleAggregation(const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                   size_t num_groups,
                   std::shared_ptr<AngularQuadrature>& quadrature,
                   std::shared_ptr<MeshContinuum>& grid);

  std::vector<AngleSetGroup> angle_set_groups;

  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& GetSimBoundaries() const
  {
    return boundaries_;
  }

  /// Resets all the outgoing intra-location and inter-location cyclic interfaces.
  void ZeroOutgoingDelayedPsi();

  /// Resets all the incoming intra-location and inter-location cyclic interfaces.
  void ZeroIncomingDelayedPsi();

  /// Initializes reflecting boundary conditions.
  void InitializeReflectingBCs();

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
};

} // namespace opensn
