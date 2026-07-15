// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/sparse_matrix/sparse_matrix.h"
#include "framework/materials/interpolator/rt_ndarray.h"
#include "framework/materials/interpolator/cartesian_grid.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include <cstdint>
#include <limits>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace opensn
{

enum XSType : std::uint64_t
{
  Total = 1 << 0,
  Absorption = 1 << 1,
  Fission = 1 << 2,
  NuFission = 1 << 3,
  Chi = 1 << 4,
  ProductionMatrix = 1 << 5,
  NuPromptFission = 1 << 6,
  NuDelayedFission = 1 << 7,
  Precursor = 1 << 8,
  Transfer = 1 << 9,
  Default = std::numeric_limits<std::uint64_t>::max(),
};

/**
 * Sparse matrix metadata.
 * Instances can be merged to form the union of the sparsity patterns found at all
 * interpolation-grid points.
 */
struct SparseMatrixProperties
{
  /// Number of matrix rows.
  std::size_t row_size;
  /// Number of matrix columns.
  std::size_t col_size;
  /// Row-col indices of each entry of the sparse matrix.
  std::set<std::pair<std::size_t, std::size_t>> present_indices;

  /// Construct properties from the dimensions and present indices of \p mat.
  SparseMatrixProperties(const SparseMatrix& mat);
  /// Add the present indices from \p prop to this sparsity pattern.
  SparseMatrixProperties& operator+=(const SparseMatrixProperties& prop);
};

/**
 * Cross-section metadata required to interpolate a MultiGroupXS.
 * This type separates invariant cross-section metadata from the values that are interpolated.
 * Transfer-matrix sparsity and delayed-neutron precursor identities may be merged across grid
 * points. These properties do not participate in equality comparison.
 */
struct XSProperties
{
  /// Number of energy groups.
  unsigned int num_groups = 0;
  /// Highest Legendre order in the transfer matrices.
  unsigned int scattering_order = 0;
  /// Number of delayed-neutron precursors in each source cross section.
  unsigned int num_precursors = 0;
  /// Flag indicating whether fission cross-section data are present.
  bool is_fissionable = false;
  /// Flag indicating whether the cross sections are in adjoint mode.
  bool adjoint_mode = false;

  /// Metadata of the transfer matrices, indexed by Legendre order.
  std::vector<SparseMatrixProperties> transfer_matrix_props;
  /// Unique precursor identities, each comprising a decay constant and delayed emission spectrum.
  std::set<std::pair<double, std::vector<double>>> precursors;
  /// Energy deposition.
  std::vector<double> energy_deposition;

  /// Construct an empty property set.
  XSProperties() = default;
  /// Extract interpolation metadata from \p xs.
  explicit XSProperties(const MultiGroupXS& xs);
  /**
   * Compare invariant metadata.
   * \note Transfer-matrix and precursor properties are not involved.
   */
  bool operator==(const XSProperties& other) const;
  /// Merge the transfer-matrix sparsity patterns in \p other into this object.
  void MergeTransferMatricesProp(const XSProperties& other);
  /// Merge precursor identities in \p other into this object, if they are not identical.
  void MergePrecursors(const XSProperties& other);

private:
  /// Insert \p precursor unless an identical identity is already present.
  void InsertPrecursor(const std::pair<double, std::vector<double>>& precursor);
};

/**
 * Base class for interpolating selected multigroup cross-section data on a Cartesian grid.
 * The constructor validates the supplied cross sections, records their common properties, and
 * packs the cross sections selected by \p flag into a contiguous array.
 *
 * Derived classes implement process coefficient computation and interpolation of the contiguous
 * array. Base class then reconstructs a MultiGroupXS from the contiguous interpolated values.
 */
class Interpolator
{
public:
  /**
   * Construct an interpolator from cross sections defined at every point of \p grid.
   *
   * \param grid Cartesian interpolation grid.
   * \param xs Cross sections vector in C-contiguous layout.
   * \param flag Bitwise combination of XSType values to interpolate.
   */
  Interpolator(const CartesianGrid& grid,
               const std::vector<std::shared_ptr<MultiGroupXS>>& xs,
               std::uint64_t flag);

  /// Return the common cross-section properties recorded during construction.
  const XSProperties& GetRefXSProp() const { return ref_xs_prop_; }
  /// Return the bitwise XSType selection used during construction.
  std::uint64_t GetFlag() const { return flag_; }

  /**
   * Evaluate the selected cross sections at \p state_point.
   * This methods check for valid in-bound state point, calls ``EvaluateContiguous``, then scatters
   * contiguous interpolated cross sections into allocated memory owned by the resulted MultiGroupXS
   * object.
   */
  MultiGroupXS Evaluate(std::span<double> state_point);

  virtual ~Interpolator() = default;

protected:
  /**
   * Get interpolated cross-sections in contiguous vector form.
   * \note \p state_point is assumed to be a valid in-bound point. The check is already performed in
   * ``Evaluate``.
   */
  virtual std::vector<double> EvaluateContiguous(std::span<double> state_point) = 0;

  XSProperties ref_xs_prop_;
  CartesianGrid grid_;
  RT_NDArray xs_data_;
  std::uint64_t flag_;
};

} // namespace opensn
