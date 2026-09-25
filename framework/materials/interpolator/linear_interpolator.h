// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/interpolator/interpolator.h"

namespace opensn
{

/**
 * Multilinear interpolation of cross-section data on a Cartesian grid.
 * For a given state point, the interpolator isolates the data at the corners of the surrounding
 * hypercube and computes the interpolation along each axis by "folding." At each fold along axis
 * \f$ x \f$, between points \f$ x_1 \f$ and \f$ x_2 \f$, the cross section \f$ \sigma \f$ is
 * computed as:
 * \f[ \sigma_1 \longleftarrow \left(1 - \frac{x-x_1}{x_2-x_1} \right) \sigma_1 +
 * \frac{x-x_1}{x_2-x_1} \sigma_2
 * \text{.} \f]
 * Each fold halves the number of cross sections. The number of folds equals the number of
 * dimensions, reducing the data at all hypercube corners to a single interpolated value. The
 * result is returned by ``EvaluateContiguous``.
 */
class LinearInterpolator : public Interpolator
{
public:
  /**
   * Construct a multilinear cross-section interpolator.
   * \param grid Cartesian interpolation grid.
   * \param xs Cross sections at every grid point, stored in C-contiguous order.
   * \param flag Bitwise combination of \ref XSType values to interpolate.
   */
  LinearInterpolator(const CartesianGrid& grid,
                     const std::vector<std::shared_ptr<MultiGroupXS>>& xs,
                     std::uint64_t flag);

  ~LinearInterpolator() override = default;

protected:
  std::vector<double> EvaluateContiguous(std::span<double> state_point) override;
};

} // namespace opensn
