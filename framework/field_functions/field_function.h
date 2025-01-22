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
public:
  static InputParameters GetInputParameters();
  explicit FieldFunction(const InputParameters& params);

  FieldFunction(std::string name, Unknown unknown);

  virtual ~FieldFunction() = default;

  /// Returns the text name of the field function.
  const std::string& Name() const { return name_; }

  /// Returns a reference to the unknown structure.
  const Unknown& GetUnknown() const { return unknown_; }

  /// Returns a reference to the unknown manager that can be with spatial discretizations.
  const UnknownManager& GetUnknownManager() const { return unknown_manager_; }

  /// Overrides the stack placement so that field functions go to the field function stack.
  void PushOntoStack(std::shared_ptr<Object> new_object) override;

  /// Evaluate the field function on a given cell, at a given position, for the given component.
  virtual double Evaluate(const Cell& cell, const Vector3& position, int component) const
  {
    return 0.0;
  }

protected:
  std::string name_;
  Unknown unknown_;
  UnknownManager unknown_manager_;
};

} // namespace opensn
