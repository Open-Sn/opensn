#pragma once

#include "framework/object.h"

namespace opensn
{

/**Base class for mesh modifiers*/
class MeshModifier : public ChiObject
{
public:
  explicit MeshModifier(const InputParameters& params);

  virtual void Apply() = 0;
  virtual ~MeshModifier() = default;

protected:
  MeshModifier() = default;
};

} // namespace opensn
