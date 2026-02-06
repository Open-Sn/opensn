// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/aah_sweep_chunk_rz.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <cmath>
#include <iomanip>
#include <stdexcept>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, DiscreteOrdinatesCurvilinearProblem);

InputParameters
DiscreteOrdinatesCurvilinearProblem::GetInputParameters()
{
  InputParameters params = DiscreteOrdinatesProblem::GetInputParameters();

  params.SetGeneralDescription(
    "Solver for Discrete Ordinates in cylindrical and spherical coordinates");

  params.SetClassName("DiscreteOrdinatesCurvilinearProblem");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesCurvilinearProblem");

  return params;
}

std::shared_ptr<DiscreteOrdinatesCurvilinearProblem>
DiscreteOrdinatesCurvilinearProblem::Create(const ParameterBlock& params)
{
  log.Log0Warning()
    << "The curvilinear discrete-ordinates problem type is experimental. USE WITH CAUTION!"
    << std::endl;
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<DiscreteOrdinatesCurvilinearProblem>(
    "lbs::DiscreteOrdinatesCurvilinearProblem", params);
}

DiscreteOrdinatesCurvilinearProblem::DiscreteOrdinatesCurvilinearProblem(
  const InputParameters& params)
  : DiscreteOrdinatesProblem(params)
{
  PerformInputChecks();
}

void
DiscreteOrdinatesCurvilinearProblem::PerformInputChecks()
{
  if (geometry_type_ != GeometryType::TWOD_CYLINDRICAL)
  {
    std::stringstream oss;
    oss << GetName() << ":\n"
        << "Invalid geometry type " << ToString(geometry_type_) << ".\n"
        << "Only TWOD_CYLINDRICAL geometry type is supported.";
    throw std::runtime_error(oss.str());
  }
  if (use_gpus_)
  {
    throw std::runtime_error(
      "DiscreteOrdinatesCurvilinearProblem: GPU acceleration is not supported for curvilinear "
      "geometries yet.");
  }

  for (size_t gs = 0; gs < groupsets_.size(); ++gs)
  {
    // angular quadrature type must be compatible with coordinate system
    const auto angular_quad_ptr = groupsets_[gs].quadrature;
    switch (grid_->GetCoordinateSystem())
    {
      case CoordinateSystemType::CYLINDRICAL:
      {
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<GLCProductQuadrature2DRZ>(angular_quad_ptr);

        if (curvilinear_angular_quad_ptr == nullptr)
        {
          std::ostringstream oss;
          oss << GetName() << ":\n"
              << "Invalid angular quadrature (type = "
              << static_cast<int>(angular_quad_ptr->GetType()) << ")";
          throw std::runtime_error(oss.str());
        }
        break;
      }
      case CoordinateSystemType::SPHERICAL:
      {
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<GLProductQuadrature1DSpherical>(angular_quad_ptr);

        if (curvilinear_angular_quad_ptr == nullptr)
        {
          std::ostringstream oss;
          oss << GetName() << ":\n"
              << "Invalid angular quadrature (type = "
              << static_cast<int>(angular_quad_ptr->GetType()) << ")";
          throw std::runtime_error(oss.str());
        }
        break;
      }
      default:
      {
        std::ostringstream oss;
        oss << GetName() << ":\n"
            << "Invalid coordinate system (type = "
            << std::to_string(static_cast<int>(grid_->GetCoordinateSystem())) << ")";
        throw std::runtime_error(oss.str());
      }
    }

    if (grid_->GetType() == MeshType::UNSTRUCTURED &&
        grid_->GetCoordinateSystem() == CoordinateSystemType::CYLINDRICAL &&
        groupsets_[gs].angleagg_method != AngleAggregationType::SINGLE)
    {
      log.Log0Warning() << GetName() << ":\n"
                        << "Forcing SINGLE angle aggregation for unstructured RZ meshes (groupset "
                        << gs << ").";
      groupsets_[gs].angleagg_method = AngleAggregationType::SINGLE;
    }

    // angle aggregation type must be compatible with coordinate system
    const auto angleagg_method = groupsets_[gs].angleagg_method;
    switch (grid_->GetCoordinateSystem())
    {
      case CoordinateSystemType::CYLINDRICAL:
      {
        if (angleagg_method != AngleAggregationType::AZIMUTHAL and
            angleagg_method != AngleAggregationType::SINGLE)
        {
          std::ostringstream oss;
          oss << GetName() << ":\n"
              << "Invalid angle aggregation (type = " << static_cast<int>(angleagg_method)
              << ") for groupsset " << gs << ". Supported: AZIMUTHAL, SINGLE.";
          throw std::runtime_error(oss.str());
        }
        break;
      }
      case CoordinateSystemType::SPHERICAL:
      {
        if (angleagg_method != AngleAggregationType::POLAR)
        {
          std::ostringstream oss;
          oss << GetName() << ":\n"
              << "Invalid angle aggregation (type = " << static_cast<int>(angleagg_method)
              << ") for groupsset " << gs;
          throw std::runtime_error(oss.str());
        }
        break;
      }
      default:
      {
        std::ostringstream oss;
        oss << GetName() << ":\nInvalid coordinate system (type = "
            << std::to_string(static_cast<int>(grid_->GetCoordinateSystem())) << ")";
        throw std::runtime_error(oss.str());
      }
    }
  }

  // boundary of mesh must be rectangular with origin at (0, 0, 0)
  const std::vector<Vector3> unit_normal_vectors = {
    Vector3(1.0, 0.0, 0.0),
    Vector3(0.0, 1.0, 0.0),
    Vector3(0.0, 0.0, 1.0),
  };
  for (const auto& cell : grid_->local_cells)
  {
    for (const auto& face : cell.faces)
    {
      if (not face.has_neighbor)
      {
        bool face_orthogonal = false;
        for (size_t d = 0; d < unit_normal_vectors.size(); ++d)
        {
          const auto n_dot_e = face.normal.Dot(unit_normal_vectors[d]);
          if (std::fabs(n_dot_e) > 0.999999)
          {
            // Allow inner radial boundaries if faces are axis-aligned
            face_orthogonal = true;
            break;
          }
        }
        if (not face_orthogonal)
        {
          std::ostringstream oss;
          oss << GetName() << ":\n"
              << "Mesh contains boundary faces not orthogonal with respect to Cartesian reference "
              << "frame";
          throw std::runtime_error(oss.str());
        }
      }
    }
  }
}

void
DiscreteOrdinatesCurvilinearProblem::InitializeSpatialDiscretization()
{
  log.Log() << "Initializing spatial discretization.\n";

  const auto quadrature_orders = [](GeometryType g) -> std::pair<QuadratureOrder, QuadratureOrder>
  {
    switch (g)
    {
      case GeometryType::ONED_SPHERICAL:
        return {QuadratureOrder::FOURTH, QuadratureOrder::THIRD};

      case GeometryType::ONED_CYLINDRICAL:
      case GeometryType::TWOD_CYLINDRICAL:
        return {QuadratureOrder::THIRD, QuadratureOrder::SECOND};

      default:
        return {QuadratureOrder::INVALID_ORDER, QuadratureOrder::INVALID_ORDER};
    }
  };

  const auto [quad_primary, quad_secondary] = quadrature_orders(geometry_type_);
  if (quad_primary == QuadratureOrder::INVALID_ORDER)
  {
    std::ostringstream oss;
    oss << GetName() << "::InitializeSpatialDiscretization:\n"
        << "Invalid geometry type " << ToString(geometry_type_) << ".\n"
        << "Only ONED_SPHERICAL, ONED_CYLINDRICAL, or TWOD_CYLINDRICAL geometries are supported.\n";
    throw std::runtime_error(oss.str());
  }

  // Primary
  discretization_ = PieceWiseLinearDiscontinuous::New(grid_, quad_primary);
  ComputeUnitIntegrals();

  // Secondary
  discretization_secondary_ = PieceWiseLinearDiscontinuous::New(grid_, quad_secondary);
  ComputeSecondaryUnitIntegrals();
}

void
DiscreteOrdinatesCurvilinearProblem::ComputeSecondaryUnitIntegrals()
{
  log.Log() << "Computing RZ secondary unit integrals.\n";
  const auto& sdm = *discretization_;

  // Secondary matrices are used for the azimuthal streaming term in RZ.
  // That term carries a 1/r factor, so use unweighted volume integrals
  // here.
  const auto swf = [](const Vector3&) { return 1.0; };

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

  const size_t num_local_cells = grid_->local_cells.size();
  secondary_unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_->local_cells)
    secondary_unit_cell_matrices_[cell.local_id] = ComputeCellUnitIntegrals(cell);

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
  auto sweep_chunk = std::make_shared<AAHSweepChunkRZ>(*this, groupset);

  return sweep_chunk;
}

} // namespace opensn
