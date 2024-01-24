#pragma once

#include "framework/math/functions/function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

/**
 * Base class for evaluating spatial material functions given a coordinate and material ID.
 */
class SpatialMaterialFunction : public Function
{
public:
  static InputParameters GetInputParameters();
  explicit SpatialMaterialFunction(const InputParameters& params);

  /**
   * Evaluate the function.
   *
   * \param xyz The xyz coordinates of the point where the function is evaluated.
   * \param mat_id Material ID
   * \param num_components Number of components
   * \return Vector with the response (should have `num_groups` entries)
   */
  virtual std::vector<double>
  Evaluate(const Vector3& xyz, int mat_id, int num_components) const = 0;
};

} // namespace opensn
