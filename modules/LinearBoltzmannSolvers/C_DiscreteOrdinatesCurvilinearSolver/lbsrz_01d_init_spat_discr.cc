#include "lbs_curvilinear_solver.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

#include <iomanip>

namespace lbs
{

void
DiscreteOrdinatesCurvilinearSolver::InitializeSpatialDiscretization()
{
  Chi::log.Log() << "Initializing spatial discretization_.\n";

  auto qorder = chi_math::QuadratureOrder::INVALID_ORDER;
  auto system = chi_math::CoordinateSystemType::UNDEFINED;

  //  primary discretisation
  switch (options_.geometry_type)
  {
    case lbs::GeometryType::ONED_SPHERICAL:
    {
      qorder = chi_math::QuadratureOrder::FOURTH;
      system = chi_math::CoordinateSystemType::SPHERICAL;
      break;
    }
    case lbs::GeometryType::ONED_CYLINDRICAL:
    case lbs::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    default:
    {
      Chi::log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::"
                                "InitializeSpatialDiscretization : "
                             << "invalid geometry, static_cast<int>(type) = "
                             << static_cast<int>(options_.geometry_type);
      Chi::Exit(EXIT_FAILURE);
    }
  }

  typedef chi_math::spatial_discretization::PieceWiseLinearDiscontinuous SDM_PWLD;
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
      qorder = chi_math::QuadratureOrder::THIRD;
      system = chi_math::CoordinateSystemType::CYLINDRICAL;
      break;
    }
    case lbs::GeometryType::ONED_CYLINDRICAL:
    case lbs::GeometryType::TWOD_CYLINDRICAL:
    {
      qorder = chi_math::QuadratureOrder::SECOND;
      system = chi_math::CoordinateSystemType::CARTESIAN;
      break;
    }
    default:
    {
      Chi::log.LogAllError() << "D_DO_RZ_SteadyState::SteadyStateSolver::"
                                "InitializeSpatialDiscretization : "
                             << "invalid geometry, static_cast<int>(type) = "
                             << static_cast<int>(options_.geometry_type);
      Chi::Exit(EXIT_FAILURE);
    }
  }

  discretization_secondary_ = SDM_PWLD::New(*grid_ptr_, qorder, system);

  ComputeSecondaryUnitIntegrals();
}

void
DiscreteOrdinatesCurvilinearSolver::ComputeSecondaryUnitIntegrals()
{
  Chi::log.Log() << "Computing RZ secondary unit integrals.\n";
  const auto& sdm = *discretization_;

  //======================================== Define spatial weighting functions
  std::function<double(const chi_mesh::Vector3&)> swf =
    chi_math::SpatialDiscretization::CylindricalRZSpatialWeightFunction;

  //======================================== Define lambda for cell-wise comps
  auto ComputeCellUnitIntegrals = [&sdm, &swf](const chi_mesh::Cell& cell)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    //    const size_t cell_num_faces = cell.faces.size();
    const size_t cell_num_nodes = cell_mapping.NumNodes();
    const auto vol_qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    MatDbl IntV_shapeI_shapeJ(cell_num_nodes, VecDbl(cell_num_nodes));

    // Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : vol_qp_data.QuadraturePointIndices())
        {
          IntV_shapeI_shapeJ[i][j] +=
            swf(vol_qp_data.QPointXYZ(qp)) * vol_qp_data.ShapeValue(i, qp) *
            vol_qp_data.ShapeValue(j, qp) * vol_qp_data.JxW(qp); // M-matrix
        }                                                        // for qp
      }                                                          // for j
    }                                                            // for i

    return lbs::UnitCellMatrices{{},                 // K-matrix
                                 {},                 // G-matrix
                                 IntV_shapeI_shapeJ, // M-matrix
                                 {},                 // Vi-vectors

                                 {},  // face M-matrices
                                 {},  // face G-matrices
                                 {}}; // face Si-vectors
  };

  const size_t num_local_cells = grid_ptr_->local_cells.size();
  secondary_unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_ptr_->local_cells)
    secondary_unit_cell_matrices_[cell.local_id_] = ComputeCellUnitIntegrals(cell);

  Chi::mpi.Barrier();
  Chi::log.Log() << "Secondary Cell matrices computed.         Process memory = "
                 << std::setprecision(3) << chi::Console::GetMemoryUsageInMB() << " MB";
}

} // namespace lbs