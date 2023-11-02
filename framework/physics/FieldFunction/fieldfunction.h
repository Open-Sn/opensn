#ifndef CHITECH_FIELDFUNCTION_H
#define CHITECH_FIELDFUNCTION_H

#include "ChiObject.h"
#include "math/UnknownManager/unknown_manager.h"
#include "mesh/chi_mesh.h"

namespace chi_mesh
{
class Cell;
}

namespace chi_physics
{

class FieldFunction : public ChiObject
{
private:
  std::string text_name_;
  chi_math::Unknown unknown_;
  chi_math::UnknownManager unknown_manager_;

public:
  /**Returns required input parameters.*/
  static chi::InputParameters GetInputParameters();

  /**ObjectMaker based constructor.*/
  explicit FieldFunction(const chi::InputParameters& params);

  /**Conventional constructor.*/
  FieldFunction(const std::string& text_name, chi_math::Unknown unknown);

  virtual ~FieldFunction() = default;

  // Getters
  /**Returns the text name of the field function.*/
  const std::string& TextName() const { return text_name_; }
  /**Returns a reference to the unknown structure.*/
  const chi_math::Unknown& Unknown() const { return unknown_; }
  /**Returns a reference to the unknown manager that can be used in
   * spatial discretizations.*/
  const chi_math::UnknownManager& GetUnknownManager() const { return unknown_manager_; }

  /**\brief Overrides the stack placement so that FieldFunctions go
   * to the field function stack.*/
  void PushOntoStack(std::shared_ptr<ChiObject>& new_object) override;

  virtual double Evaluate(const chi_mesh::Cell& cell,
                          const chi_mesh::Vector3& position,
                          unsigned int component) const
  {
    return 0.0;
  }
};

} // namespace chi_physics

#endif // CHITECH_FIELDFUNCTION_H
