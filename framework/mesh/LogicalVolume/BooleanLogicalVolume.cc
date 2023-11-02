#include "BooleanLogicalVolume.h"

#include "ChiObjectFactory.h"

namespace chi_mesh
{

chi::InputParameters BooleanLogicalVolumeArgumentPair();

RegisterChiObject(chi_mesh, BooleanLogicalVolume);
RegisterSyntaxBlock(chi_mesh, BooleanLogicalVolumeArgumentPair, BooleanLogicalVolumeArgumentPair);

chi::InputParameters
BooleanLogicalVolume::GetInputParameters()
{
  chi::InputParameters params = LogicalVolume::GetInputParameters();

  // clang-format off
  params.SetDocGroup("LuaLogicVolumes\n");
  // clang-format on

  params.AddRequiredParameterArray(
    "parts",
    "Array of combinatorial logic each entry has the following required params "
    "<TT>chi_mesh::BooleanLogicalVolumeArgumentPair</TT>");

  params.LinkParameterToBlock("parts", "chi_mesh::BooleanLogicalVolumeArgumentPair");

  return params;
}

BooleanLogicalVolume::BooleanLogicalVolume(const chi::InputParameters& params)
  : LogicalVolume(params)
{
  const auto& input_parts = params.GetParam("parts");
  input_parts.RequireBlockTypeIs(chi::ParameterBlockType::ARRAY);

  for (size_t p = 0; p < input_parts.NumParameters(); ++p)
  {
    const auto& part = input_parts.GetParam(p);
    part.RequireBlockTypeIs(chi::ParameterBlockType::BLOCK);

    auto part_params = BooleanLogicalVolumeArgumentPair();

    part_params.AssignParameters(part);

    const size_t lv_handle = part_params.GetParamValue<size_t>("lv");
    auto lv_ptr =
      Chi::GetStackItemPtrAsType<LogicalVolume>(Chi::object_stack, lv_handle, __FUNCTION__);

    parts.emplace_back(part_params.GetParamValue<bool>("op"), lv_ptr);
  }
}

chi::InputParameters
BooleanLogicalVolumeArgumentPair()
{
  chi::InputParameters params;

  params.SetDocGroup("chi_mesh__BooleanLogicalVolume");

  params.AddRequiredParameter<bool>(
    "op",
    "Boolean value indicating the volume sense. True means inside, False means "
    "outside");
  params.AddRequiredParameter<size_t>("lv", "Handle to a logical volume.");

  return params;
}

bool
BooleanLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  for (const auto& part : parts)
  {
    if (part.first != part.second->Inside(point)) return false;
  }

  return true;
}

} // namespace chi_mesh