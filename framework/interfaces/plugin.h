#pragma once

#include "framework/object.h"

namespace opensn
{

class Plugin : public ChiObject
{
public:
  static InputParameters GetInputParameters();
  explicit Plugin(const InputParameters& params);

  ~Plugin();

protected:
  const std::string plugin_path_;
  void* library_handle_ = nullptr;
};

} // namespace opensn
