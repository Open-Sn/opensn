// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/mat.h"
#include "framework/logging/log.h"
#include <memory>

using namespace opensn;

namespace opensnlua
{

std::shared_ptr<opensn::Material>
MatAddMaterial(const std::string& name)
{
  auto new_material = std::make_shared<Material>();
  new_material->name = name;

  opensn::material_stack.push_back(new_material);

  const size_t index = opensn::material_stack.size() - 1;

  opensn::log.Log0Verbose1() << "New material added at index " << index << " with name \""
                             << new_material->name << "\"";
  return new_material;
}

} // namespace opensnlua
