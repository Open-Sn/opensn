#include "lua_mesh_modifiers.h"

#include "framework/runtime.h"

#include "framework/mesh/mesh_modifiers/mesh_modifier.h"

#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterWrapperFunction(chi_mesh,
                        MeshModifiersApply,
                        MeshModifiersApply_Syntax,
                        MeshModifiersApply);

InputParameters
MeshModifiersApply_Syntax()
{
  InputParameters params;

  params.SetGeneralDescription("Lua wrapper function for applying mesh modifiers");
  params.SetDocGroup("DocMeshModifiers");

  params.AddRequiredParameterArray("arg0", "A list of handles to the modifiers to apply.");

  return params;
}

ParameterBlock
MeshModifiersApply(const InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const std::vector<size_t> handles = params.GetParamVectorValue<size_t>("arg0");

  for (const size_t handle : handles)
  {
    auto& modifier =
      opensn::Chi::GetStackItem<MeshModifier>(opensn::Chi::object_stack, handle, fname);

    modifier.Apply();
  }

  return ParameterBlock(); // Return empty param block
}

} // namespace opensnlua
