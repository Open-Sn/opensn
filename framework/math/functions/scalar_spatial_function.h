#pragma once

#include "framework/math/functions/function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

/**
 * Base class for evaluating functions given spatial location
 */
class ScalarSpatialFunction : public Function
{
public:
  static InputParameters GetInputParameters();
  explicit ScalarSpatialFunction(const InputParameters& params);

  /**
   * Evaluate this function
   *
   * \param xyz The xyz coordinates of the point where the function is called.
   * \return Function value
   */
  virtual double Evaluate(const Vector3& xyz) const = 0;
};

} // namespace opensn
