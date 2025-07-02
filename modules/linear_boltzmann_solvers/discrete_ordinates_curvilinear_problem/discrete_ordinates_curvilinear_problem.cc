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
#include <iomanip>

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
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesCurvilinearProblem");

  return params;
}

std::shared_ptr<DiscreteOrdinatesCurvilinearProblem>
DiscreteOrdinatesCurvilinearProblem::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<DiscreteOrdinatesCurvilinearProblem>(
    "lbs::DiscreteOrdinatesCurvilinearProblem", params);
}

DiscreteOrdinatesCurvilinearProblem::DiscreteOrdinatesCurvilinearProblem(
  const InputParameters& params)
  : DiscreteOrdinatesProblem(params)
{
}

void
DiscreteOrdinatesCurvilinearProblem::PerformInputChecks()
{
  log.Log() << "DiscreteOrdinatesCurvilinearSolver::PerformInputChecks : enter";

  DiscreteOrdinatesProblem::PerformInputChecks();

  // perform additional verifications for curvilinear LBS
  // coordinate system must be curvilinear
  if (grid_->GetCoordinateSystem() != CoordinateSystemType::CYLINDRICAL and
      grid_->GetCoordinateSystem() != CoordinateSystemType::SPHERICAL)
  {
    std::ostringstream oss;
    oss << "DiscreteOrdinatesCurvilinearSolver: Invalid coordinate system (type = "
        << std::to_string(static_cast<int>(grid_->GetCoordinateSystem())) << ")";
    throw std::runtime_error(oss.str());
  }

  // re-interpret geometry type to curvilinear
  switch (options_.geometry_type)
  {
    case GeometryType::ONED_SLAB:
    {
      std::ostringstream oss;
      oss << "DiscreteOrdinatesCurvilinearProblem: Invalid geometry (type = "
          << std::to_string(static_cast<int>(grid_->GetCoordinateSystem())) << ")";
      throw std::runtime_error(oss.str());
      break;
    }
    case GeometryType::TWOD_CARTESIAN:
    {
      switch (grid_->GetCoordinateSystem())
      {
        case CoordinateSystemType::CYLINDRICAL:
        {
          options_.geometry_type = GeometryType::TWOD_CYLINDRICAL;
          break;
        }
        default:
        {
          std::ostringstream oss;
          oss << "DiscreteOrdinatesCurvilinearSolver: Invalid geometry (type = "
              << std::to_string(static_cast<int>(options_.geometry_type)) << ") "
              << "for curvilinear coordinate system (type = "
              << std::to_string(static_cast<int>(grid_->GetCoordinateSystem())) << ")";
          throw std::runtime_error(oss.str());
        }
      }
      break;
    }
    default:
    {
      std::ostringstream oss;
      oss << "DiscreteOrdinatesCurvilinearProblem: Invalid geometry (type = "
          << std::to_string(static_cast<int>(options_.geometry_type)) << ") "
          << "for curvilinear coordinate system";
      throw std::runtime_error(oss.str());
    }
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
          oss << "DiscreteOrdinatesCurvilinearProblem: Invalid angular quadrature (type = "
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
          oss << "DiscreteOrdinatesCurvilinearProblem: Invalid angular quadrature (type = "
              << static_cast<int>(angular_quad_ptr->GetType()) << ")";
          throw std::runtime_error(oss.str());
        }
        break;
      }
      default:
      {
        std::ostringstream oss;
        oss << "DiscreteOrdinatesCurvilinearProblem: Invalid coordinate system (type = "
            << std::to_string(static_cast<int>(grid_->GetCoordinateSystem())) << ")";
        throw std::runtime_error(oss.str());
      }
    }

    // angle aggregation type must be compatible with coordinate system
    const auto angleagg_method = groupsets_[gs].angleagg_method;
    switch (grid_->GetCoordinateSystem())
    {
      case CoordinateSystemType::CYLINDRICAL:
      {
        if (angleagg_method != AngleAggregationType::AZIMUTHAL)
        {
          std::ostringstream oss;
          oss << "DiscreteOrdinatesCurvilinearProblem: Invalid angle aggregation (type = "
              << static_cast<int>(angleagg_method) << ") for groupsset " << gs;
          throw std::runtime_error(oss.str());
        }
        break;
      }
      case CoordinateSystemType::SPHERICAL:
      {
        if (angleagg_method != AngleAggregationType::POLAR)
        {
          std::ostringstream oss;
          oss << "DiscreteOrdinatesCurvilinearProblem: Invalid angle aggregation (type = "
              << static_cast<int>(angleagg_method) << ") for groupsset " << gs;
          throw std::runtime_error(oss.str());
        }
        break;
      }
      default:
      {
        std::ostringstream oss;
        oss << "DiscreteOrdinatesCurvilinearProblem: Invalid coordinate system (type = "
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
          if (n_dot_e > 0.999999)
          {
            face_orthogonal = true;
            break;
          }
          else if (n_dot_e < -0.999999)
          {
            for (const auto& v_id : face.vertex_ids)
            {
              const auto& vertex = grid_->vertices[v_id];
              if (std::abs(vertex[d]) > 1.0e-12)
              {
                std::ostringstream oss;
                oss << "DiscreteOrdinatesCurvilinearProblem: Mesh contains boundary faces with "
                    << "outward-oriented unit normal vector " +
                         (-1 * unit_normal_vectors[d]).PrintStr()
                    << ", with vertices characterized by v(" + std::to_string(d) + ") != 0";
                throw std::runtime_error(oss.str());
              }
            }
            face_orthogonal = true;
            break;
          }
        }
        if (not face_orthogonal)
        {
          std::ostringstream oss;
          oss << "DiscreteOrdinatesCurvilinearProblem: Mesh contains boundary faces not "
              << "orthogonal with respect to Cartesian reference frame";
          throw std::runtime_error(oss.str());
        }
      }
    }
  }
  log.Log() << "DiscreteOrdinatesCurvilinearSolver::PerformInputChecks : exit";
}

void
DiscreteOrdinatesCurvilinearProblem::InitializeSpatialDiscretization()
{
  log.Log() << "Initializing spatial discretization_.\n";

  // primary discretization
  QuadratureOrder qorder;
  switch (options_.geometry_type)
  {
    case GeometryType::ONED_SPHERICAL:
    {
      qorder = QuadratureOrder::FOURTH;
      break;
    }
    case GeometryType::ONED_CYLINDRICAL:
    case GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = QuadratureOrder::THIRD;
      break;
    }
    default:
    {
      std::ostringstream oss;
      oss << "DiscreteOrdinatesCurvilinearProblem: Invalid geometry (type = "
          << static_cast<int>(options_.geometry_type) << ")";
      throw std::runtime_error(oss.str());
    }
  }

  discretization_ = PieceWiseLinearDiscontinuous::New(grid_, qorder);

  ComputeUnitIntegrals();

  // secondary discretization
  // system - manipulated such that the spatial discretization returns
  // a cell view of the same type but with weighting of degree one less
  // than the primary discretization
  switch (options_.geometry_type)
  {
    case GeometryType::ONED_SPHERICAL:
    {
      qorder = QuadratureOrder::THIRD;
      break;
    }
    case GeometryType::ONED_CYLINDRICAL:
    case GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = QuadratureOrder::SECOND;
      break;
    }
    default:
    {
      std::ostringstream oss;
      oss << "DiscreteOrdinatesCurvilinearProblem: Invalid geometry (type = "
          << static_cast<int>(options_.geometry_type) << ")";
      throw std::runtime_error(oss.str());
    }
  }

  discretization_secondary_ = PieceWiseLinearDiscontinuous::New(grid_, qorder);

  ComputeSecondaryUnitIntegrals();
}

void
DiscreteOrdinatesCurvilinearProblem::ComputeSecondaryUnitIntegrals()
{
  log.Log() << "Computing RZ secondary unit integrals.\n";
  const auto& sdm = *discretization_;

  // Define spatial weighting functions
  std::function<double(const Vector3&)> swf =
    SpatialDiscretization::CylindricalRZSpatialWeightFunction;

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
      }   // for j
    }     // for i

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

std::shared_ptr<SweepChunk>
DiscreteOrdinatesCurvilinearProblem::SetSweepChunk(LBSGroupset& groupset)
{
  auto sweep_chunk = std::make_shared<AAHSweepChunkRZ>(grid_,
                                                       *discretization_,
                                                       unit_cell_matrices_,
                                                       secondary_unit_cell_matrices_,
                                                       cell_transport_views_,
                                                       densities_local_,
                                                       phi_new_local_[groupset.id],
                                                       psi_new_local_[groupset.id],
                                                       q_moments_local_[groupset.id],
                                                       groupset,
                                                       block_id_to_xs_map_,
                                                       num_moments_,
                                                       max_cell_dof_count_);

  return sweep_chunk;
}

} // namespace opensn
