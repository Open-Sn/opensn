#include "framework/mesh/FieldFunctionInterpolation/Point/chi_ffinter_point.h"

#include "framework/physics/FieldFunction/fieldfunction_gridbased.h"
#include "framework/math/SpatialDiscretization/SpatialDiscretization.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/chi_runtime.h"
#include "framework/mpi/chi_mpi.h"

// ###################################################################
/**Initializes the point interpolator.*/
void
chi_mesh::FieldFunctionInterpolationPoint::Initialize()
{
  const std::string fname = "FieldFunctionInterpolationPoint::Initialize";
  const auto& grid = field_functions_.front()->GetSpatialDiscretization().Grid();

  std::vector<uint64_t> cells_potentially_owning_point;
  for (const auto& cell : grid.local_cells)
  {
    const auto& vcc = cell.centroid_;
    const auto& poi = point_of_interest_;
    const auto nudged_point = poi + 1.0e-6 * (vcc - poi);
    if (grid.CheckPointInsideCell(cell, nudged_point))
      cells_potentially_owning_point.push_back(cell.global_id_);
  }

  const int local_count = static_cast<int>(cells_potentially_owning_point.size());
  std::vector<int> locI_count(Chi::mpi.process_count, 0);
  MPI_Allgather(&local_count, // sendbuf
                1,
                MPI_INT,           // sendcount + sendtype
                locI_count.data(), // recvbuf + recvtype
                1,
                MPI_INT,        // recvcount + recvtype
                Chi::mpi.comm); // communicator

  std::vector<int> recvdispls(Chi::mpi.process_count, 0);

  int running_count = 0;
  for (int locI = 0; locI < Chi::mpi.process_count; ++locI)
  {
    recvdispls[locI] = running_count;
    running_count += locI_count[locI];
  }

  const auto& sendbuf = cells_potentially_owning_point;
  std::vector<uint64_t> recvbuf(running_count, 0);
  MPI_Allgatherv(sendbuf.data(), // sendbuf
                 local_count,
                 MPI_UINT64_T,      // sendcount + sendtype
                 recvbuf.data(),    // recvbuf
                 locI_count.data(), // recvcount
                 recvdispls.data(), // recvdispl
                 MPI_UINT64_T,      // recvtype
                 Chi::mpi.comm);    // communicator

  if (recvbuf.empty()) throw std::logic_error(fname + ": No cell identified containing the point.");

  uint64_t owning_cell_gid = recvbuf.front();
  for (const uint64_t gid : recvbuf)
    owning_cell_gid = std::min(owning_cell_gid, gid);

  locally_owned_ = false;
  for (const uint64_t gid : cells_potentially_owning_point)
    if (gid == owning_cell_gid)
    {
      locally_owned_ = true;
      owning_cell_gid_ = owning_cell_gid;
      break;
    }
}

// ###################################################################
/**Executes the point interpolator.*/
void
chi_mesh::FieldFunctionInterpolationPoint::Execute()
{
  if (not locally_owned_) return;

  const auto& ref_ff = *field_functions_.front();
  const auto& sdm = ref_ff.GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.GetUnknownManager();
  const auto uid = 0;
  const auto cid = ref_component_;

  using namespace chi_mesh::ff_interpolation;
  const auto field_data = ref_ff.GetGhostedFieldVector();

  const auto& cell = grid.cells[owning_cell_gid_];
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t num_nodes = cell_mapping.NumNodes();

  std::vector<double> node_dof_values(num_nodes, 0.0);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
    node_dof_values[i] = field_data[imap];
  } // for i

  std::vector<double> shape_values(num_nodes, 0.0);
  cell_mapping.ShapeValues(point_of_interest_, shape_values);

  point_value_ = 0.0;
  for (size_t i = 0; i < num_nodes; ++i)
    point_value_ += node_dof_values[i] * shape_values[i];
}

// ###################################################################
/**Gets the value of the field function evaluation at the point.*/
double
chi_mesh::FieldFunctionInterpolationPoint::GetPointValue() const
{
  double global_point_value;
  MPI_Allreduce(&point_value_,       // sendbuf
                &global_point_value, // recvbuf
                1,
                MPI_DOUBLE,     // count + datatype
                MPI_SUM,        // operation
                Chi::mpi.comm); // communicator

  return global_point_value;
}
