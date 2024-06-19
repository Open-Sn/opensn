// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/mesh/mesh.h"

namespace opensn
{
class Cell;

class FieldFunction : public Object
{
private:
  std::string text_name_;
  Unknown unknown_;
  UnknownManager unknown_manager_;

public:
  /**Returns required input parameters.*/
  static InputParameters GetInputParameters();

  /**ObjectMaker based constructor.*/
  explicit FieldFunction(const InputParameters& params);

  /**Conventional constructor.*/
  FieldFunction(const std::string& text_name, Unknown unknown);

  virtual ~FieldFunction() = default;

  // Getters
  /**Returns the text name of the field function.*/
  const std::string& TextName() const { return text_name_; }
  /**Returns a reference to the unknown structure.*/
  const Unknown& GetUnknown() const { return unknown_; }
  /**Returns a reference to the unknown manager that can be used in
   * spatial discretizations.*/
  const UnknownManager& GetUnknownManager() const { return unknown_manager_; }

  /**\brief Overrides the stack placement so that FieldFunctions go
   * to the field function stack.*/
  void PushOntoStack(std::shared_ptr<Object>& new_object) override;

  virtual double Evaluate(const Cell& cell, const Vector3& position, unsigned int component) const
  {
    return 0.0;
  }
};

} // namespace opensn
