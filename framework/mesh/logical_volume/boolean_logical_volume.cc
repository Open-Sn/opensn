#include "framework/mesh/logical_volume/boolean_logical_volume.h"

#include "framework/object_factory.h"

namespace opensn
{

InputParameters BooleanLogicalVolumeArgumentPair();

OpenSnRegisterObjectInNamespace(mesh, BooleanLogicalVolume);
OpenSnRegisterSyntaxBlockInNamespace(mesh,
                                     BooleanLogicalVolumeArgumentPair,
                                     BooleanLogicalVolumeArgumentPair);

InputParameters
BooleanLogicalVolume::GetInputParameters()
{
  InputParameters params = LogicalVolume::GetInputParameters();

  // clang-format off
  params.SetDocGroup("LuaLogicVolumes\n");
  // clang-format on

  params.AddRequiredParameterArray(
    "parts",
    "Array of combinatorial logic each entry has the following required params "
    "<TT>mesh::BooleanLogicalVolumeArgumentPair</TT>");

  params.LinkParameterToBlock("parts", "mesh::BooleanLogicalVolumeArgumentPair");

  return params;
}

BooleanLogicalVolume::BooleanLogicalVolume(const InputParameters& params) : LogicalVolume(params)
{
  const auto& input_parts = params.GetParam("parts");
  input_parts.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  for (size_t p = 0; p < input_parts.NumParameters(); ++p)
  {
    const auto& part = input_parts.GetParam(p);
    part.RequireBlockTypeIs(ParameterBlockType::BLOCK);

    auto part_params = BooleanLogicalVolumeArgumentPair();

    part_params.AssignParameters(part);

    const size_t lv_handle = part_params.GetParamValue<size_t>("lv");
    auto lv_ptr = GetStackItemPtrAsType<LogicalVolume>(object_stack, lv_handle, __FUNCTION__);

    parts.emplace_back(part_params.GetParamValue<bool>("op"), lv_ptr);
  }
}

InputParameters
BooleanLogicalVolumeArgumentPair()
{
  InputParameters params;

  params.SetDocGroup("mesh__BooleanLogicalVolume");

  params.AddRequiredParameter<bool>(
    "op",
    "Boolean value indicating the volume sense. True means inside, False means "
    "outside");
  params.AddRequiredParameter<size_t>("lv", "Handle to a logical volume.");

  return params;
}

bool
BooleanLogicalVolume::Inside(const Vector3& point) const
{
  for (const auto& part : parts)
  {
    if (part.first != part.second->Inside(point))
      return false;
  }

  return true;
}

} // namespace opensn
