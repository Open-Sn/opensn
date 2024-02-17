#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_curvilinear_solver/lbs_curvilinear_solver.h"
#include "framework/object_factory.h"
#include "framework/math/quadratures/cylindrical_angular_quadrature.h"
#include "framework/math/quadratures/spherical_angular_quadrature.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_curvilinear_solver/sweep_chunks/lbs_curvilinear_sweep_chunk_pwl.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/memory_usage.h"
#include <iomanip>

namespace opensn
{
namespace lbs
{

OpenSnRegisterObjectInNamespace(lbs, DiscreteOrdinatesCurvilinearSolver);

InputParameters
DiscreteOrdinatesCurvilinearSolver::GetInputParameters()
{
  InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();

  params.SetGeneralDescription(
    "Solver for Discrete Ordinates in cylindrical and spherical coordinates");

  params.SetClassName("DiscreteOrdinatesCurvilinearSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesCurvilinearSolver");
  params.AddRequiredParameter<int>("coord_system",
                                   "Coordinate system to use. Can only be 2 or "
                                   "3. 2=Cylindrical, 3=Spherical.");

  return params;
}

DiscreteOrdinatesCurvilinearSolver::DiscreteOrdinatesCurvilinearSolver(
  const InputParameters& params)
  : DiscreteOrdinatesSolver(params),
    coord_system_type_(static_cast<CoordinateSystemType>(params.GetParamValue<int>("coord_system")))
{
}

void
DiscreteOrdinatesCurvilinearSolver::PerformInputChecks()
{
  log.Log() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : enter";

  //  --------------------------------------------------------------------------
  //  perform all verifications of Cartesian LBS
  //  --------------------------------------------------------------------------

  lbs::DiscreteOrdinatesSolver::PerformInputChecks();

  //  --------------------------------------------------------------------------
  //  perform additional verifications for curvilinear LBS
  //  --------------------------------------------------------------------------

  //  coordinate system must be curvilinear
  if (coord_system_type_ != CoordinateSystemType::CYLINDRICAL and
      coord_system_type_ != CoordinateSystemType::SPHERICAL)
  {
    log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                      << "invalid coordinate system, static_cast<int>(type) = "
                      << static_cast<int>(coord_system_type_);
    Exit(EXIT_FAILURE);
  }

  //  re-interpret geometry type to curvilinear
  switch (options_.geometry_type)
  {
    case lbs::GeometryType::ONED_SLAB:
    {
      switch (coord_system_type_)
      {
        default:
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "invalid geometry, static_cast<int>(type) = "
                            << static_cast<int>(options_.geometry_type) << " "
                            << "for curvilinear coordinate system, static_cast<int>(type) = "
                            << static_cast<int>(coord_system_type_);
          Exit(EXIT_FAILURE);
        }
      }
      break;
    }
    case lbs::GeometryType::TWOD_CARTESIAN:
    {
      switch (coord_system_type_)
      {
        case CoordinateSystemType::CYLINDRICAL:
        {
          options_.geometry_type = lbs::GeometryType::TWOD_CYLINDRICAL;
          break;
        }
        default:
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "invalid geometry, static_cast<int>(type) = "
                            << static_cast<int>(options_.geometry_type) << " "
                            << "for curvilinear coordinate system, static_cast<int>(type) = "
                            << static_cast<int>(coord_system_type_);
          Exit(EXIT_FAILURE);
        }
      }
      break;
    }
    default:
    {
      log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                        << "invalid geometry, static_cast<int>(type) = "
                        << static_cast<int>(options_.geometry_type) << " "
                        << "for curvilinear coordinate system";
      Exit(EXIT_FAILURE);
    }
  }

  for (size_t gs = 0; gs < groupsets_.size(); ++gs)
  {
    //  angular quadrature type must be compatible with coordinate system
    const auto angular_quad_ptr = groupsets_[gs].quadrature_;
    switch (coord_system_type_)
    {
      case CoordinateSystemType::CYLINDRICAL:
      {
        typedef CylindricalAngularQuadrature CylAngQuad;
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<CylAngQuad>(angular_quad_ptr);
        if (curvilinear_angular_quad_ptr == nullptr)
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "invalid angular quadrature, static_cast<int>(type) = "
                            << static_cast<int>(angular_quad_ptr->type_)
                            << ", for groupset = " << gs;
          Exit(EXIT_FAILURE);
        }
        break;
      }
      case CoordinateSystemType::SPHERICAL:
      {
        typedef SphericalAngularQuadrature SphAngQuad;
        const auto curvilinear_angular_quad_ptr =
          std::dynamic_pointer_cast<SphAngQuad>(angular_quad_ptr);
        if (curvilinear_angular_quad_ptr == nullptr)
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "invalid angular quadrature, static_cast<int>(type) = "
                            << static_cast<int>(angular_quad_ptr->type_)
                            << ", for groupset = " << gs;
          Exit(EXIT_FAILURE);
        }
        break;
      }
      default:
      {
        log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
                          << static_cast<int>(coord_system_type_);
        Exit(EXIT_FAILURE);
      }
    }

    //  angle aggregation type must be compatible with coordinate system
    const auto angleagg_method = groupsets_[gs].angleagg_method_;
    switch (coord_system_type_)
    {
      case CoordinateSystemType::CYLINDRICAL:
      {
        if (angleagg_method != lbs::AngleAggregationType::AZIMUTHAL)
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "invalid angle aggregation type, static_cast<int>(type) = "
                            << static_cast<int>(angleagg_method) << ", for groupset = " << gs;
          Exit(EXIT_FAILURE);
        }
        break;
      }
      case CoordinateSystemType::SPHERICAL:
      {
        if (angleagg_method != lbs::AngleAggregationType::POLAR)
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "invalid angle aggregation type, static_cast<int>(type) = "
                            << static_cast<int>(angleagg_method) << ", for groupset = " << gs;
          Exit(EXIT_FAILURE);
        }
        break;
      }
      default:
      {
        log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                          << "invalid curvilinear coordinate system, static_cast<int>(type) = "
                          << static_cast<int>(coord_system_type_);
        Exit(EXIT_FAILURE);
      }
    }
  }

  //  boundary of mesh must be rectangular with origin at (0, 0, 0)
  const std::vector<Vector3> unit_normal_vectors = {
    Vector3(1.0, 0.0, 0.0),
    Vector3(0.0, 1.0, 0.0),
    Vector3(0.0, 0.0, 1.0),
  };
  for (const auto& cell : grid_ptr_->local_cells)
  {
    for (const auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
      {
        bool face_orthogonal = false;
        for (size_t d = 0; d < unit_normal_vectors.size(); ++d)
        {
          const auto n_dot_e = face.normal_.Dot(unit_normal_vectors[d]);
          if (n_dot_e > 0.999999)
          {
            face_orthogonal = true;
            break;
          }
          else if (n_dot_e < -0.999999)
          {
            for (const auto& v_id : face.vertex_ids_)
            {
              const auto& vertex = grid_ptr_->vertices[v_id];
              if (std::abs(vertex[d]) > 1.0e-12)
              {
                log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::"
                                     "PerformInputChecks : "
                                  << "mesh contains boundary faces with outward-oriented unit "
                                  << "normal vector " << (-1 * unit_normal_vectors[d]).PrintS()
                                  << "with vertices characterised by v(" << d << ") != 0.";
                Exit(EXIT_FAILURE);
              }
            }
            face_orthogonal = true;
            break;
          }
        }
        if (not face_orthogonal)
        {
          log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : "
                            << "mesh contains boundary faces not orthogonal with respect to "
                            << "Cartesian reference frame.";
          Exit(EXIT_FAILURE);
        }
      }
    }
  }

  log.Log() << "D_DO_RZ_SteadyState::SteadyStateSolver::PerformInputChecks : exit";
}

void
DiscreteOrdinatesCurvilinearSolver::InitializeSpatialDiscretization()
{
  log.Log() << "Initializing spatial discretization_.\n";

  auto qorder = QuadratureOrder::INVALID_ORDER;
  auto system = CoordinateSystemType::UNDEFINED;

  //  primary discretisation
  switch (options_.geometry_type)
  {
    case lbs::GeometryType::ONED_SPHERICAL:
    {
      qorder = QuadratureOrder::FOURTH;
      system = CoordinateSystemType::SPHERICAL;
      break;
    }
    case lbs::GeometryType::ONED_CYLINDRICAL:
    case lbs::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = QuadratureOrder::THIRD;
      system = CoordinateSystemType::CYLINDRICAL;
      break;
    }
    default:
    {
      log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::"
                           "InitializeSpatialDiscretization : "
                        << "invalid geometry, static_cast<int>(type) = "
                        << static_cast<int>(options_.geometry_type);
      Exit(EXIT_FAILURE);
    }
  }

  typedef PieceWiseLinearDiscontinuous SDM_PWLD;
  discretization_ = SDM_PWLD::New(*grid_ptr_, qorder, system);

  ComputeUnitIntegrals();

  //  secondary discretisation
  //  system - manipulated such that the spatial discretisation returns
  //  a cell view of the same type but with weighting of degree one less
  //  than the primary discretisation
  switch (options_.geometry_type)
  {
    case lbs::GeometryType::ONED_SPHERICAL:
    {
      qorder = QuadratureOrder::THIRD;
      system = CoordinateSystemType::CYLINDRICAL;
      break;
    }
    case lbs::GeometryType::ONED_CYLINDRICAL:
    case lbs::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = QuadratureOrder::SECOND;
      system = CoordinateSystemType::CARTESIAN;
      break;
    }
    default:
    {
      log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::"
                           "InitializeSpatialDiscretization : "
                        << "invalid geometry, static_cast<int>(type) = "
                        << static_cast<int>(options_.geometry_type);
      Exit(EXIT_FAILURE);
    }
  }

  discretization_secondary_ = SDM_PWLD::New(*grid_ptr_, qorder, system);

  ComputeSecondaryUnitIntegrals();
}

void
DiscreteOrdinatesCurvilinearSolver::ComputeSecondaryUnitIntegrals()
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
    const size_t cell_num_nodes = cell_mapping.NumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    MatDbl IntV_shapeI_shapeJ(cell_num_nodes, VecDbl(cell_num_nodes));

    // Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : fe_vol_data.QuadraturePointIndices())
        {
          IntV_shapeI_shapeJ[i][j] += swf(fe_vol_data.QPointXYZ(qp)) *
                                      fe_vol_data.ShapeValue(i, qp) *
                                      fe_vol_data.ShapeValue(j, qp) * fe_vol_data.JxW(qp);
        } // for qp
      }   // for j
    }     // for i

    return lbs::UnitCellMatrices{{},
                                 {},
                                 IntV_shapeI_shapeJ,
                                 {},

                                 {},
                                 {},
                                 {}};
  };

  const size_t num_local_cells = grid_ptr_->local_cells.size();
  secondary_unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_ptr_->local_cells)
    secondary_unit_cell_matrices_[cell.local_id_] = ComputeCellUnitIntegrals(cell);

  opensn::mpi_comm.barrier();
  log.Log() << "Secondary Cell matrices computed.         Process memory = " << std::setprecision(3)
            << GetMemoryUsageInMB() << " MB";
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesCurvilinearSolver::SetSweepChunk(lbs::LBSGroupset& groupset)
{
  auto sweep_chunk = std::make_shared<SweepChunkPwlrz>(*grid_ptr_,
                                                       *discretization_,
                                                       unit_cell_matrices_,
                                                       secondary_unit_cell_matrices_,
                                                       cell_transport_views_,
                                                       phi_new_local_,
                                                       psi_new_local_[groupset.id_],
                                                       q_moments_local_,
                                                       groupset,
                                                       matid_to_xs_map_,
                                                       num_moments_,
                                                       max_cell_dof_count_);

  return sweep_chunk;
}

} // namespace lbs
} // namespace opensn
