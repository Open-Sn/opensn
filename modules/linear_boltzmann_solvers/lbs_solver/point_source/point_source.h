// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh.h"
#include "framework/object.h"
#include <utility>

namespace opensn
{

class LBSSolver;

/**
 * A class for point sources, which is defined by its location and a group-wise strength vector.
 *
 * A point source can belong to one or more cells based on whether it lies on the interior of a
 * cell, a face, or a vertex. When a point source lies on a face or vertex, it belongs to multiple
 * cells. Each of these cells are called subscribers. In this case, the source strength is split
 * amongst the cells using a volumetric weighting. On each cell, the point source is defined in the
 * finite element basis to obtain the contribution to each node. Defining the point source as
 * \f$ Q(x) = S \delta(x - x_0) = \sum_j s_j b_j \f$, this implies that \f[ \int Q(x) dx = S \int
 * \delta(x - x_0) dx = \sum_j s_j \int b_j dx. \f] To solve for the coefficients, a Galerkin
 * method is used such that \f[ S \int b_i \delta(x - x_0) dx = \sum_j s_j \int b_i b_j dx. \f]
 * Using the standard definition of the mass matrix, this is given by \f[ S \vec{b}(x_0) = M
 * \vec{s}. \f] This can be solved for the coefficients via
 * \f[ \vec{s} = S M^{-1} \vec{b}(x_0). \f].
 */
class PointSource : public Object
{
public:
  /// An encapsulation of subscriber cell data for defining point source contributions.
  struct Subscriber
  {
    double volume_weight;
    uint64_t cell_local_id;
    std::vector<double> shape_values;
    std::vector<double> node_weights;
  };

public:
  static InputParameters GetInputParameters();
  explicit PointSource(const InputParameters& params);

  /**
   * Initializes the cell subscriber info from the given solver.
   *
   * @note Uninitialized point sources will not contribute to sources.
   */
  void Initialize(const LBSSolver& lbs_solver);

  size_t NumLocalSubscribers() const { return subscribers_.size(); }
  size_t NumGlobalSubscribers() const { return num_global_subscribers_; }
  const std::vector<Subscriber>& Subscribers() const { return subscribers_; }

  const Vector3& Location() const { return location_; }
  const std::vector<double>& Strength() const { return strength_; }

private:
  const Vector3 location_;

  // TODO: This could be a map with group indices to avoid
  //       potentially many zero entries.
  const std::vector<double> strength_;

  std::vector<Subscriber> subscribers_;
  size_t num_global_subscribers_ = 0;
};

} // namespace opensn
