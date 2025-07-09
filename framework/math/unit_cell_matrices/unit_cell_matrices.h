// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/dense_matrix.h"
#include "framework/math/geometry.h"
#include "framework/mesh/cell/cell.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/vector3.h"

namespace opensn
{

struct UnitCellMatrices
{
  DenseMatrix<double> intV_gradshapeI_gradshapeJ;
  DenseMatrix<Vector3> intV_shapeI_gradshapeJ;
  DenseMatrix<double> intV_shapeI_shapeJ;
  Vector<double> intV_shapeI;

  std::vector<DenseMatrix<double>> intS_shapeI_shapeJ;
  std::vector<DenseMatrix<Vector3>> intS_shapeI_gradshapeJ;
  std::vector<Vector<double>> intS_shapeI;

  /// Compute unit cell matrices for a given cell
  static UnitCellMatrices
  Compute(const SpatialDiscretization& sdm, const Cell& cell, GeometryType geometry_type);

  /// Compute unit cell matrices for a given cell assuming cartesian geometry
  static UnitCellMatrices Compute(const SpatialDiscretization& sdm, const Cell& cell);
};

} // namespace opensn
