// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/logical_volume/boolean_logical_volume.h"
#include "framework/mesh/mesh.h"
#include "framework/object_factory.h"

namespace opensn
{

InputParameters BooleanLogicalVolumeArgumentPair();

OpenSnRegisterObjectInNamespace(logvol, BooleanLogicalVolume);

InputParameters
BooleanLogicalVolume::GetInputParameters()
{
  InputParameters params = LogicalVolume::GetInputParameters();

  params.AddRequiredParameterArray(
    "parts",
    "Array of combinatorial logic each entry has the following required params "
    "<TT>mesh::BooleanLogicalVolumeArgumentPair</TT>");

  params.LinkParameterToBlock("parts", "mesh::BooleanLogicalVolumeArgumentPair");

  return params;
}

std::shared_ptr<BooleanLogicalVolume>
BooleanLogicalVolume::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<BooleanLogicalVolume>("logvol::BooleanLogicalVolume", params);
}

BooleanLogicalVolume::BooleanLogicalVolume(const InputParameters& params) : LogicalVolume(params)
{
  const auto& input_parts = params.GetParam("parts");
  input_parts.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  for (size_t p = 0; p < input_parts.GetNumParameters(); ++p)
  {
    const auto& part = input_parts.GetParam(p);
    part.RequireBlockTypeIs(ParameterBlockType::BLOCK);

    auto part_params = BooleanLogicalVolumeArgumentPair();

    part_params.AssignParameters(part);

    auto lv = part_params.GetSharedPtrParam<LogicalVolume>("lv", false);
    parts.emplace_back(part_params.GetParamValue<bool>("op"), lv);
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
  params.AddRequiredParameter<std::shared_ptr<LogicalVolume>>("lv", "Logical volume.");

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
