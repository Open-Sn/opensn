#pragma once

#include "framework/math/functions/function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

/**
 * Base class for evaluating response functions given spatial location and material ID
 */
class ResponseFunction : public Function
{
public:
  static InputParameters GetInputParameters();
  explicit ResponseFunction(const InputParameters& params);

  /**
   * Evaluate this function
   *
   * \param num_groups Number of groups
   * \param xyz The xyz coordinates of the point where the function is evaluated.
   * \param mat_id Material ID
   * \return Vector with the response (should have `num_groups` entries)
   */
  virtual std::vector<double> Evaluate(int num_groups, const Vector3& xyz, int mat_id) const = 0;
};

} // namespace opensn
