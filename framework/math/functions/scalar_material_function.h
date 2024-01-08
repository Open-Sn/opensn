#pragma once

#include "framework/math/functions/function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

/**
 * Base class for evaluating functions given material ID and single value
 */
class ScalarMaterialFunction : public Function
{
public:
  static InputParameters GetInputParameters();
  explicit ScalarMaterialFunction(const InputParameters& params);

  /**
   * Evaluate this function
   *
   * \param val The scalar value (for example, a field function value)
   * \param mat_id The material ID of the cell
   * \return Function value
   */
  virtual double Evaluate(double val, int mat_id) const = 0;
};

} // namespace opensn
