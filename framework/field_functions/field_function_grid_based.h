// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/field_function.h"

#include "framework/mesh/mesh.h"

#include "framework/math/parallel_vector/ghosted_parallel_stl_vector.h"
#include <string>
#include <memory>
#include <vector>
#include <utility>

namespace opensn
{
class SpatialDiscretization;
class GhostedParallelSTLVector;

/***/
class FieldFunctionGridBased : public FieldFunction
{
public:
  typedef std::pair<Vector3, Vector3> BoundingBox;

  /**Returns required input parameters.*/
  static InputParameters GetInputParameters();

  /**ObjectMaker based constructor.*/
  explicit FieldFunctionGridBased(const InputParameters& params);

  /**Creates a field function, filling it with zeros.*/
  FieldFunctionGridBased(const std::string& text_name,
                         std::shared_ptr<SpatialDiscretization>& discretization_ptr,
                         Unknown unknown);

  /**Creates a field function with an associated field vector.
   * The field's data vector is set to the incoming field vector.*/
  FieldFunctionGridBased(const std::string& text_name,
                         std::shared_ptr<SpatialDiscretization>& sdm_ptr,
                         Unknown unknown,
                         const std::vector<double>& field_vector);

  /**Creates a field function where all the values are assigned to
   * the single supplied value.*/
  FieldFunctionGridBased(const std::string& text_name,
                         std::shared_ptr<SpatialDiscretization>& sdm_ptr,
                         Unknown unknown,
                         double field_value);

  virtual ~FieldFunctionGridBased() = default;

  /**
   * Returns the spatial discretization method.
   */
  const SpatialDiscretization& GetSpatialDiscretization() const;

  /**
   * Returns a read-only reference to the locally stored field data.
   */
  const std::vector<double>& FieldVectorRead() const;

  /**
   * Returns a reference to the locally stored field data.
   */
  std::vector<double>& FieldVector();

  /**
   * Updates the field vector with a local STL vector.
   */
  void UpdateFieldVector(const std::vector<double>& field_vector);

  /**
   * Updates the field vector with a PETSc vector. This only operates locally.
   */
  void UpdateFieldVector(const Vec& field_vector);

  /**
   * Static method to export multiple grid-based field functions.
   */
  typedef std::vector<std::shared_ptr<const FieldFunctionGridBased>> FFList;

  /**
   * Export multiple field functions to VTK.
   */
  static void ExportMultipleToVTK(const std::string& file_base_name, const FFList& ff_list);

  /**
   * Makes a copy of the locally stored data with ghost access.
   */
  std::vector<double> GetGhostedFieldVector() const;

  /**
   * Returns the component values at requested point.
   */
  virtual std::vector<double> GetPointValue(const Vector3& point) const;

  /**Evaluates the field function, on a cell, at the specified point.*/
  double Evaluate(const Cell& cell, const Vector3& position, unsigned int component) const override;

protected:
  std::shared_ptr<SpatialDiscretization> sdm_;
  std::unique_ptr<GhostedParallelSTLVector> ghosted_field_vector_;

private:
  /**
   * Static method for making the GetSpatialDiscretization for the constructors.
   */
  static std::shared_ptr<SpatialDiscretization>
  MakeSpatialDiscretization(const InputParameters& params);

  /**
   * Static method for making the ghosted vector for the constructors.
   */
  static std::unique_ptr<GhostedParallelSTLVector>

  /**
   * Private method for creating the field vector.
   */
  MakeFieldVector(const SpatialDiscretization& discretization, const UnknownManager& uk_man);

  const BoundingBox local_grid_bounding_box_;
};

} // namespace opensn
