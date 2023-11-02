#include "lua_meshmodifiers.h"

#include "chi_runtime.h"

#include "mesh/MeshModifiers/MeshModifier.h"

#include "console/chi_console.h"

namespace chi_mesh::lua_utils
{

RegisterWrapperFunction(/*namespace_in_lua=*/chi_mesh,
                        /*name_in_lua=*/MeshModifiersApply,
                        /*syntax_function=*/MeshModifiersApply_Syntax,
                        /*actual_function=*/MeshModifiersApply);

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