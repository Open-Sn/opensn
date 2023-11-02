#include "SphereLogicalVolume.h"

#include "ChiObjectFactory.h"

namespace chi_mesh
{

RegisterChiObject(chi_mesh, SphereLogicalVolume);

chi::InputParameters
SphereLogicalVolume::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetDocGroup("LuaLogicVolumes");

  params.AddOptionalParameter("r", 1.0, "Radius of the sphere.");
  params.AddOptionalParameter("x", 0.0, "X-location of the volume.");
  params.AddOptionalParameter("y", 0.0, "Y-location of the volume.");
  params.AddOptionalParameter("z", 0.0, "Z-location of the volume.");

  using namespace chi_data_types;
  params.ConstrainParameterRange("r", AllowableRangeLowLimit::New(0.0, /*low_closed=*/false));

  return params;
}

SphereLogicalVolume::SphereLogicalVolume(const chi::InputParameters& params)
  : LogicalVolume(params),
    r_(params.GetParamValue<double>("r")),
    x0_(params.GetParamValue<double>("x")),
    y0_(params.GetParamValue<double>("y")),
    z0_(params.GetParamValue<double>("z"))
{
}

bool
SphereLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  double dx = point.x - x0_;
  double dy = point.y - y0_;
  double dz = point.z - z0_;

  double R2 = dx * dx + dy * dy + dz * dz;

  if (R2 <= (r_ * r_)) return true;
  else
    return false;
}

} // namespace chi_mesh