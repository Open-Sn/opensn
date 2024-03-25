#include "framework/mesh/logical_volume/sphere_logical_volume.h"

#include "framework/object_factory.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(logvol, SphereLogicalVolume);

InputParameters
SphereLogicalVolume::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetDocGroup("LuaLogicVolumes");

  params.AddOptionalParameter("r", 1.0, "Radius of the sphere.");
  params.AddOptionalParameter("x", 0.0, "X-location of the volume.");
  params.AddOptionalParameter("y", 0.0, "Y-location of the volume.");
  params.AddOptionalParameter("z", 0.0, "Z-location of the volume.");

  params.ConstrainParameterRange("r", AllowableRangeLowLimit::New(0.0, false));

  return params;
}

SphereLogicalVolume::SphereLogicalVolume(const InputParameters& params)
  : LogicalVolume(params),
    r_(params.GetParamValue<double>("r")),
    x0_(params.GetParamValue<double>("x")),
    y0_(params.GetParamValue<double>("y")),
    z0_(params.GetParamValue<double>("z"))
{
}

bool
SphereLogicalVolume::Inside(const Vector3& point) const
{
  double dx = point.x - x0_;
  double dy = point.y - y0_;
  double dz = point.z - z0_;

  double R2 = dx * dx + dy * dy + dz * dz;

  if (R2 <= (r_ * r_))
    return true;
  else
    return false;
}

} // namespace opensn
