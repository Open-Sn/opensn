#include "fieldfunction_gridbased.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"

namespace chi_physics
{

// ##################################################################
/**Not necessarily an optimal routine to use when repeatedly querying
 * a field function.*/
std::vector<double>
FieldFunctionGridBased::GetPointValue(const chi_mesh::Vector3& point) const
{
  typedef const int64_t cint64_t;
  const auto& uk_man = GetUnknownManager();
  const size_t num_components = uk_man.GetTotalUnknownStructureSize();

  size_t local_num_point_hits = 0;
  std::vector<double> local_point_value(num_components, 0.0);

  const auto& xyz_min = local_grid_bounding_box_.first;
  const auto& xyz_max = local_grid_bounding_box_.second;

  const double xmin = xyz_min.x;
  const double ymin = xyz_min.y;
  const double zmin = xyz_min.z;

  const double xmax = xyz_max.x;
  const double ymax = xyz_max.y;
  const double zmax = xyz_max.z;

  const auto& field_vector = *ghosted_field_vector_;

  if (point.x >= xmin and point.x <= xmax and point.y >= ymin and point.y <= ymax and
      point.z >= zmin and point.z <= zmax)
  {
    const auto& grid = sdm_->Grid();
    for (const auto& cell : grid.local_cells)
    {
      if (grid.CheckPointInsideCell(cell, point))
      {
        const auto& cell_mapping = sdm_->GetCellMapping(cell);
        std::vector<double> shape_values;
        cell_mapping.ShapeValues(point, shape_values);

        local_num_point_hits += 1;

        const size_t num_nodes = cell_mapping.NumNodes();
        for (size_t c = 0; c < num_components; ++c)
        {
          for (size_t j = 0; j < num_nodes; ++j)
          {
            cint64_t dof_map_j = sdm_->MapDOFLocal(cell, j, uk_man, 0, c);
            const double dof_value_j = field_vector[dof_map_j];

            local_point_value[c] += dof_value_j * shape_values[j];
          } // for node i
        }   // for component c
      }     // if inside cell
    }       // for cell
  }         // if in bounding box

  //============================================= Communicate number of
  //                                              point hits
  size_t globl_num_point_hits;
  MPI_Allreduce(&local_num_point_hits, // sendbuf
                &globl_num_point_hits, // recvbuf
                1,
                MPIU_SIZE_T,    // count + datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  std::vector<double> globl_point_value(num_components, 0.0);
  MPI_Allreduce(local_point_value.data(), // sendbuf
                globl_point_value.data(), // recvbuf
                1,
                MPI_DOUBLE,     // count + datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  chi_math::Scale(globl_point_value, 1.0 / static_cast<double>(globl_num_point_hits));

  return globl_point_value;
}

// ##################################################################
double
FieldFunctionGridBased::Evaluate(const chi_mesh::Cell& cell,
                                 const chi_mesh::Vector3& position,
                                 unsigned int component) const
{
  const auto& field_vector = *ghosted_field_vector_;

  typedef const int64_t cint64_t;
  const auto& cell_mapping = sdm_->GetCellMapping(cell);

  std::vector<double> shape_values;
  cell_mapping.ShapeValues(position, shape_values);

  double value = 0.0;
  const size_t num_nodes = cell_mapping.NumNodes();
  for (size_t j = 0; j < num_nodes; ++j)
  {
    cint64_t dof_map = sdm_->MapDOFLocal(cell, j, GetUnknownManager(), 0, component);

    value += field_vector[dof_map] * shape_values[j];
  }

  return value;
}

} // namespace chi_physics