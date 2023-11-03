#pragma once

#include "ChiObject.h"

namespace chi_mesh
{

/**Base class for mesh modifiers*/
class MeshModifier : public ChiObject
{
public:
  explicit MeshModifier(const chi::InputParameters& params);

  virtual void Apply() = 0;
  virtual ~MeshModifier() = default;

protected:
  MeshModifier() = default;
};

} // namespace chi_mesh
