#include "lua_meshmodifiers.h"

#include "framework/runtime.h"

#include "framework/mesh/mesh_modifiers/mesh_modifier.h"

#include "framework/console/console.h"

namespace chi_mesh::lua_utils
{

RegisterWrapperFunction(chi_mesh,
                        MeshModifiersApply,
                        MeshModifiersApply_Syntax,
                        MeshModifiersApply);

chi::InputParameters
MeshModifiersApply_Syntax()
{
  chi::InputParameters params;

  params.SetGeneralDescription("Lua wrapper function for applying mesh modifiers");
  params.SetDocGroup("DocMeshModifiers");

  params.AddRequiredParameterArray("arg0", "A list of handles to the modifiers to apply.");

  return params;
}

chi::ParameterBlock
MeshModifiersApply(const chi::InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const std::vector<size_t> handles = params.GetParamVectorValue<size_t>("arg0");

  for (const size_t handle : handles)
  {
    auto& modifier = Chi::GetStackItem<MeshModifier>(Chi::object_stack, handle, fname);

    modifier.Apply();
  }

  return chi::ParameterBlock(); // Return empty param block
}

} // namespace chi_mesh::lua_utils
