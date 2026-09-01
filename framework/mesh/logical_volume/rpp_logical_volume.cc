// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/parameters/input_parameters.h"
#include "framework/utils/error.h"

namespace opensn
{

InputParameters
RPPLogicalVolume::GetInputParameters()
{
  InputParameters params = LogicalVolume::GetInputParameters();

  params.AddOptionalParameter("xmin", 0.0, "X-min of the volume");
  params.AddOptionalParameter("xmax", 1.0, "X-max of the volume");
  params.AddOptionalParameter("ymin", 0.0, "Y-min of the volume");
  params.AddOptionalParameter("ymax", 1.0, "Y-max of the volume");
  params.AddOptionalParameter("zmin", 0.0, "Z-min of the volume");
  params.AddOptionalParameter("zmax", 1.0, "Z-max of the volume");

  params.AddOptionalParameter("infx", false, "Flag, when true, will ignore xmin and xmax.");
  params.AddOptionalParameter("infy", false, "Flag, when true, will ignore ymin and ymax.");
  params.AddOptionalParameter("infz", false, "Flag, when true, will ignore zmin and zmax.");

  return params;
}

std::shared_ptr<RPPLogicalVolume>
RPPLogicalVolume::Create(const ParameterBlock& params)
{
  return CreateObject<RPPLogicalVolume>("logvol::RPPLogicalVolume", params);
}

RPPLogicalVolume::RPPLogicalVolume(const InputParameters& params)
  : LogicalVolume(params),
    xmin_(params.GetParamValue<double>("xmin")),
    xmax_(params.GetParamValue<double>("xmax")),
    ymin_(params.GetParamValue<double>("ymin")),
    ymax_(params.GetParamValue<double>("ymax")),
    zmin_(params.GetParamValue<double>("zmin")),
    zmax_(params.GetParamValue<double>("zmax")),
    infx_(params.GetParamValue<bool>("infx")),
    infy_(params.GetParamValue<bool>("infy")),
    infz_(params.GetParamValue<bool>("infz"))
{
  // An inverted range (e.g. xmin > xmax) silently defines an empty volume rather than the
  // swapped range a typo would suggest -- reject it.
  OpenSnInvalidArgumentIf(not infx_ and xmin_ > xmax_, "xmin must not be greater than xmax.");
  OpenSnInvalidArgumentIf(not infy_ and ymin_ > ymax_, "ymin must not be greater than ymax.");
  OpenSnInvalidArgumentIf(not infz_ and zmin_ > zmax_, "zmin must not be greater than zmax.");
}

bool
RPPLogicalVolume::Inside(const Vector3& point) const
{
  constexpr std::array<bool, 6> true_condition = {true, true, true, true, true, true};
  std::array<bool, 6> condition = {false, false, false, false, false, false};

  if (point.x <= xmax_ or infx_)
    condition[XMAX] = true;
  if (point.x >= xmin_ or infx_)
    condition[XMIN] = true;
  if (point.y <= ymax_ or infy_)
    condition[YMAX] = true;
  if (point.y >= ymin_ or infy_)
    condition[YMIN] = true;
  if (point.z <= zmax_ or infz_)
    condition[ZMAX] = true;
  if (point.z >= zmin_ or infz_)
    condition[ZMIN] = true;

  return condition == true_condition;
}

} // namespace opensn
