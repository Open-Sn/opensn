#include "fieldfunction_gridbased.h"

#include "math/SpatialDiscretization/FiniteVolume/FiniteVolume.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "math/SpatialDiscretization/FiniteElement/Lagrange/LagrangeContinuous.h"
#include "math/SpatialDiscretization/FiniteElement/Lagrange/LagrangeDiscontinuous.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"

namespace chi_physics
{

RegisterChiObject(chi_physics, FieldFunctionGridBased);

chi::InputParameters FieldFunctionGridBased::GetInputParameters()
{
  chi::InputParameters params = FieldFunction::GetInputParameters();

  params.SetDocGroup("DocFieldFunction");

  params.AddRequiredParameter<std::string>(
    "sdm_type", "The spatial discretization type to be used");

  params.AddOptionalParameter(
    "initial_value", 0.0, "The initial value to assign to the field function");

  params.AddOptionalParameter(
    "quadrature_order",
    0,
    "If supplied, will overwrite the default for the "
    "specific discretization-coordinate system combination.");

  params.AddOptionalParameter("coordinate_system",
                              "cartesian",
                              "Coordinate system to apply to element mappings");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "sdm_type",
    AllowableRangeList::New({"FV", "PWLC", "PWLD", "LagrangeC", "LagrangeD"}));
  params.ConstrainParameterRange(
    "coordinate_system",
    AllowableRangeList::New({"cartesian", "cylindrical", "spherical"}));

  return params;
}

// ##################################################################
/**ObjectMaker based constructor.*/
FieldFunctionGridBased::FieldFunctionGridBased(
  const chi::InputParameters& params)
  : FieldFunction(params),
    sdm_(MakeSpatialDiscretization(params)),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(
      chi_mesh::GetCurrentHandler().GetGrid()->GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(params.GetParamValue<double>("initial_value"));
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(
  const std::string& text_name,
  chi_math::SDMPtr& discretization_ptr,
  chi_math::Unknown unknown)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(discretization_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(
  const std::string& text_name,
  chi_math::SDMPtr& sdm_ptr,
  chi_math::Unknown unknown,
  const std::vector<double>& field_vector)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  ChiInvalidArgumentIf(
    field_vector.size() != ghosted_field_vector_->LocalSize(),
    "Constructor called with incompatible size field vector.");

  ghosted_field_vector_->Set(field_vector);
}

// ##################################################################
FieldFunctionGridBased::FieldFunctionGridBased(const std::string& text_name,
                                               chi_math::SDMPtr& sdm_ptr,
                                               chi_math::Unknown unknown,
                                               double field_value)
  : FieldFunction(text_name, std::move(unknown)),
    sdm_(sdm_ptr),
    ghosted_field_vector_(MakeFieldVector(*sdm_, GetUnknownManager())),
    local_grid_bounding_box_(sdm_->Grid().GetLocalBoundingBox())
{
  ghosted_field_vector_->Set(field_value);
}

// ##################################################################
/**Returns the spatial discretization method.*/
const chi_math::SpatialDiscretization&
FieldFunctionGridBased::GetSpatialDiscretization() const
{
  return *sdm_;
}

// ##################################################################
/**Returns a read-only reference to the locally stored field data.*/
const std::vector<double>& FieldFunctionGridBased::FieldVectorRead() const
{
  return ghosted_field_vector_->LocalSTLData();
}
/**Returns a reference to the locally stored field data.*/
std::vector<double>& FieldFunctionGridBased::FieldVector()
{
  return ghosted_field_vector_->LocalSTLData();
}

// ##################################################################
/**Private method for creating the spatial discretization method.*/
chi_math::SDMPtr FieldFunctionGridBased::MakeSpatialDiscretization(
  const chi::InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();
  const auto& grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
  const auto sdm_type = params.GetParamValue<std::string>("sdm_type");

  typedef chi_math::spatial_discretization::FiniteVolume FV;
  typedef chi_math::spatial_discretization::PieceWiseLinearContinuous PWLC;
  typedef chi_math::spatial_discretization::PieceWiseLinearDiscontinuous PWLD;
  typedef chi_math::spatial_discretization::LagrangeContinuous LagC;
  typedef chi_math::spatial_discretization::LagrangeDiscontinuous LagD;

  if (sdm_type == "FV") return FV::New(*grid_ptr);

  chi_math::CoordinateSystemType cs_type =
    chi_math::CoordinateSystemType::CARTESIAN;
  std::string cs = "cartesian";
  if (user_params.Has("coordinate_system"))
  {
    cs = params.GetParamValue<std::string>("coordinate_system");

    using namespace chi_math;
    if (cs == "cartesian") cs_type = CoordinateSystemType::CARTESIAN;
    if (cs == "cylindrical") cs_type = CoordinateSystemType::CYLINDRICAL;
    if (cs == "spherical") cs_type = CoordinateSystemType::SPHERICAL;
  }

  chi_math::QuadratureOrder q_order = chi_math::QuadratureOrder::SECOND;

  if (user_params.Has("quadrature_order"))
  {
    const uint32_t max_order =
      static_cast<uint32_t>(chi_math::QuadratureOrder::FORTYTHIRD);
    const uint32_t q_order_int =
      params.GetParamValue<uint32_t>("quadrature_order");
    ChiInvalidArgumentIf(q_order_int > max_order,
                         "Invalid quadrature point order " +
                           std::to_string(q_order_int));
    q_order = static_cast<chi_math::QuadratureOrder>(q_order_int);
  }
  else // Defaulted
  {
    if (cs == "cartesian") q_order = chi_math::QuadratureOrder::SECOND;
    if (cs == "cylindrical") q_order = chi_math::QuadratureOrder::THIRD;
    if (cs == "spherical") q_order = chi_math::QuadratureOrder::FOURTH;
  }

  if (sdm_type == "PWLC") return PWLC::New(*grid_ptr, q_order, cs_type);
  else if (sdm_type == "PWLD")
    return PWLD::New(*grid_ptr, q_order, cs_type);
  else if (sdm_type == "LagrangeC")
    return LagC::New(*grid_ptr, q_order, cs_type);
  else if (sdm_type == "LagrangeD")
    return LagD::New(*grid_ptr, q_order, cs_type);

  // If not returned by now
  ChiInvalidArgument("Unsupported sdm_type \"" + sdm_type + "\"");
}

/**Private method for creating the field vector.*/
std::unique_ptr<chi_math::GhostedParallelSTLVector>
FieldFunctionGridBased::MakeFieldVector(
  const chi_math::SpatialDiscretization& discretization,
  const chi_math::UnknownManager& uk_man)
{
  auto field = std::make_unique<chi_math::GhostedParallelSTLVector>(
    discretization.GetNumLocalDOFs(uk_man),
    discretization.GetNumGlobalDOFs(uk_man),
    discretization.GetGhostDOFIndices(uk_man),
    Chi::mpi.comm);

  return field;
}

} // namespace chi_physics