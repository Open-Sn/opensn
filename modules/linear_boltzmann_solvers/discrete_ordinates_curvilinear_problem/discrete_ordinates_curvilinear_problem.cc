// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/aah_sweep_chunk_rz.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/aah_sweep_chunk_spherical.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_definition.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace opensn
{

InputParameters
DiscreteOrdinatesCurvilinearProblem::GetInputParameters()
{
  InputParameters params = DiscreteOrdinatesProblem::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesCurvilinearProblem");

  return params;
}

std::shared_ptr<DiscreteOrdinatesCurvilinearProblem>
DiscreteOrdinatesCurvilinearProblem::Create(const ParameterBlock& params)
{
  return Build(
    std::shared_ptr<DiscreteOrdinatesCurvilinearProblem>(new DiscreteOrdinatesCurvilinearProblem(
      MakeInputParameters<DiscreteOrdinatesCurvilinearProblem>(
        "lbs::DiscreteOrdinatesCurvilinearProblem", params))));
}

DiscreteOrdinatesCurvilinearProblem::DiscreteOrdinatesCurvilinearProblem(
  const InputParameters& params)
  : DiscreteOrdinatesProblem(params)
{
}

void
DiscreteOrdinatesCurvilinearProblem::CheckConfigurationErrors(
  std::vector<std::string>& errors) const
{
  DiscreteOrdinatesProblem::CheckConfigurationErrors(errors);

  if (geometry_type_ != GeometryType::TWOD_CYLINDRICAL and
      geometry_type_ != GeometryType::ONED_SPHERICAL)
    errors.emplace_back("Invalid geometry type " + std::string(ToString(geometry_type_)) +
                        ". Only TWOD_CYLINDRICAL and ONED_SPHERICAL geometry types are supported.");
  if (sweep_type_ != "AAH")
    errors.emplace_back("Curvilinear geometries support only sweep_type=\"AAH\".");
  if (use_gpus_)
    errors.emplace_back("GPU acceleration is not supported for curvilinear geometries yet.");

  const auto coordinate_system = grid_->GetCoordinateSystem();
  const bool cylindrical = coordinate_system == CoordinateSystemType::CYLINDRICAL;
  const bool spherical = coordinate_system == CoordinateSystemType::SPHERICAL;
  if (not cylindrical and not spherical)
    errors.emplace_back("Invalid coordinate system (type = " +
                        std::to_string(static_cast<int>(coordinate_system)) + ").");

  for (const auto& groupset : groupsets_)
  {
    const auto gs = std::to_string(groupset.id);
    // The angular quadrature must match the coordinate system.
    const auto& quadrature = groupset.quadrature;
    if (quadrature and
        ((cylindrical and not std::dynamic_pointer_cast<GLCProductQuadrature2DRZ>(quadrature)) or
         (spherical and not std::dynamic_pointer_cast<GLProductQuadrature1DSpherical>(quadrature))))
      errors.emplace_back("Invalid angular quadrature (type = " +
                          std::to_string(static_cast<int>(quadrature->GetType())) +
                          ") for groupset " + gs + ".");

    // The angle aggregation must match the coordinate system.
    const auto aggregation = groupset.angleagg_method;
    if (cylindrical and aggregation != AngleAggregationType::AZIMUTHAL and
        aggregation != AngleAggregationType::SINGLE)
      errors.emplace_back(
        "Invalid angle aggregation (type = " + std::to_string(static_cast<int>(aggregation)) +
        ") for groupset " + gs + ". Supported: AZIMUTHAL, SINGLE.");
    if (cylindrical and grid_->GetType() == MeshType::UNSTRUCTURED and
        aggregation != AngleAggregationType::SINGLE)
      errors.emplace_back("Groupset " + gs +
                          ": unstructured RZ meshes require angle_aggregation_type \"single\".");
    if (spherical and aggregation != AngleAggregationType::AZIMUTHAL and
        aggregation != AngleAggregationType::SINGLE)
      errors.emplace_back("Groupset " + gs +
                          ": 1D spherical geometry requires angle_aggregation_type \"azimuthal\" "
                          "(one inward and one outward angle set) or \"single\".");
  }

  // Every boundary face must be orthogonal to a Cartesian axis. Collective.
  const std::vector<Vector3> unit_normal_vectors = {
    Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(0.0, 0.0, 1.0)};
  int local_non_orthogonal = 0;
  for (const auto& cell : grid_->GetLocalCells())
    for (const auto& face : cell->faces)
      if (not face.has_neighbor and
          std::none_of(unit_normal_vectors.begin(),
                       unit_normal_vectors.end(),
                       [&face](const Vector3& e)
                       { return std::fabs(face.normal.Dot(e)) > 0.999999; }))
        local_non_orthogonal = 1;
  int non_orthogonal = 0;
  mpi_comm.all_reduce(local_non_orthogonal, non_orthogonal, mpi::op::max<int>());
  if (non_orthogonal != 0)
    errors.emplace_back(
      "Mesh contains boundary faces not orthogonal with respect to Cartesian reference frame.");

  // A reflecting outer surface is valid in 1D spherical geometry (mu -> -mu); in RZ the outer
  // radial boundary has no mirror direction in the quadrature.
  const auto& boundary_names = grid_->GetBoundaryNameMap();
  const auto rmax = boundary_names.find("xmax");
  if (cylindrical and rmax != boundary_names.end())
  {
    const auto definition = boundary_definitions_.find(rmax->second);
    if (definition != boundary_definitions_.end() and
        definition->second.type == LBSBoundaryType::REFLECTING)
      errors.emplace_back("Reflecting boundary on rmax is not supported in RZ. Please use vacuum "
                          "or isotropic on rmax.");
  }
}

void
DiscreteOrdinatesCurvilinearProblem::InitializeSpatialDiscretization()
{
  log.Log() << "Initializing spatial discretization.\n";

  // The radial weight raises the polynomial order of the integrands. The secondary integrals (one
  // power of r less) use the same quadrature, which also integrates them exactly.
  auto quad_order = QuadratureOrder::INVALID_ORDER;
  switch (geometry_type_)
  {
    case GeometryType::ONED_SPHERICAL:
      quad_order = QuadratureOrder::FOURTH;
      break;
    case GeometryType::ONED_CYLINDRICAL:
    case GeometryType::TWOD_CYLINDRICAL:
      quad_order = QuadratureOrder::THIRD;
      break;
    default:
      break;
  }
  // CheckConfigurationErrors has rejected other geometries before the build.
  OpenSnLogicalErrorIf(quad_order == QuadratureOrder::INVALID_ORDER,
                       GetName() + ": unsupported geometry " +
                         std::string(ToString(geometry_type_)) + ".");

  discretization_ = PieceWiseLinearDiscontinuous::New(grid_, quad_order);
  ComputeUnitIntegrals();
  ComputeSecondaryUnitIntegrals();
}

void
DiscreteOrdinatesCurvilinearProblem::ComputeSecondaryUnitIntegrals()
{
  log.Log() << "Computing curvilinear secondary unit integrals.\n";
  const auto& sdm = *discretization_;

  // Secondary matrices are used for the angular-derivative terms, which carry a 1/r factor. The
  // primary integrals use the radial weight (r in RZ, r^2 in 1D spherical geometry), so these use
  // one power of r less: 1 in RZ and r = z in 1D spherical geometry.
  const bool spherical = grid_->GetCoordinateSystem() == CoordinateSystemType::SPHERICAL;
  const auto swf = [spherical](const Vector3& pt) { return spherical ? pt.z : 1.0; };

  // Define lambda for cell-wise comps
  auto ComputeCellUnitIntegrals = [&sdm, &swf](const Cell& cell)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    //    const size_t cell_num_faces = cell.faces.size();
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    DenseMatrix<double> IntV_shapeI_shapeJ(cell_num_nodes, cell_num_nodes, 0.0);

    // Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
        {
          IntV_shapeI_shapeJ(i, j) += swf(fe_vol_data.QPointXYZ(qp)) *
                                      fe_vol_data.ShapeValue(i, qp) *
                                      fe_vol_data.ShapeValue(j, qp) * fe_vol_data.JxW(qp);
        } // for qp
      } // for j
    } // for i

    return UnitCellMatrices{{},
                            {},
                            IntV_shapeI_shapeJ,
                            {},

                            {},
                            {},
                            {}};
  };

  const size_t num_local_cells = grid_->GetLocalCellCount();
  secondary_unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_->GetLocalCells())
    secondary_unit_cell_matrices_[cell->local_id] = ComputeCellUnitIntegrals(*cell);

  opensn::mpi_comm.barrier();
  log.Log() << "Secondary Cell matrices computed.";
}

const std::vector<UnitCellMatrices>&
DiscreteOrdinatesCurvilinearProblem::GetSecondaryUnitCellMatrices() const
{
  return secondary_unit_cell_matrices_;
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesCurvilinearProblem::SetSweepChunk(LBSGroupset& groupset)
{
  if (geometry_type_ == GeometryType::ONED_SPHERICAL)
    return std::make_shared<AAHSweepChunkSpherical>(*this, groupset);
  return std::make_shared<AAHSweepChunkRZ>(*this, groupset);
}

} // namespace opensn
