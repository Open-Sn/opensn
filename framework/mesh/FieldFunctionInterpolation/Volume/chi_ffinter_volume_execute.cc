#include "chi_ffinter_volume.h"

#include "math/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Executes the volume interpolation.*/
void
chi_mesh::FieldFunctionInterpolationVolume::Execute()
{
  const auto& ref_ff = *field_functions_.front();
  const auto& sdm = ref_ff.GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.GetUnknownManager();
  const auto uid = 0;
  const auto cid = ref_component_;

  using namespace chi_mesh::ff_interpolation;
  const auto field_data = ref_ff.GetGhostedFieldVector();

  double local_volume = 0.0;
  double local_sum = 0.0;
  double local_max = 0.0;
  double local_min = 0.0;
  for (const uint64_t cell_local_id : cell_local_ids_inside_logvol_)
  {
    const auto& cell = grid.local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    std::vector<double> node_dof_values(num_nodes, 0.0);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
      node_dof_values[i] = field_data[imap];
    } // for i

    if (cell_local_id == cell_local_ids_inside_logvol_.front())
    {
      local_max = node_dof_values.front();
      local_min = node_dof_values.front();
    }

    for (size_t i = 0; i < num_nodes; ++i)
    {
      local_max = std::fmax(node_dof_values[i], local_max);
      local_min = std::fmin(node_dof_values[i], local_min);
    }

    for (const size_t qp : qp_data.QuadraturePointIndices())
    {
      double ff_value = 0.0;
      for (size_t j = 0; j < num_nodes; ++j)
        ff_value += qp_data.ShapeValue(j, qp) * node_dof_values[j];

      double function_value = ff_value;
      if (op_type_ >= Operation::OP_SUM_LUA and op_type_ <= Operation::OP_MAX_LUA)
        function_value = CallLuaFunction(ff_value, cell.material_id_);

      local_volume += qp_data.JxW(qp);
      local_sum += function_value * qp_data.JxW(qp);
      local_max = std::fmax(ff_value, local_max);
      local_min = std::fmin(ff_value, local_min);
    } // for qp
  }   // for cell-id

  if (op_type_ == Operation::OP_SUM or op_type_ == Operation::OP_SUM_LUA)
  {
    double global_sum;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, Chi::mpi.comm);
    op_value_ = global_sum;
  }
  if (op_type_ == Operation::OP_AVG or op_type_ == Operation::OP_AVG_LUA)
  {
    double local_data[] = {local_volume, local_sum};
    double global_data[] = {0.0, 0.0};

    MPI_Allreduce(&local_data, &global_data, 2, MPI_DOUBLE, MPI_SUM, Chi::mpi.comm);
    double global_volume = global_data[0];
    double global_sum = global_data[1];
    op_value_ = global_sum / global_volume;
  }
  if (op_type_ == Operation::OP_MAX or op_type_ == Operation::OP_MAX_LUA)
  {
    double global_value;
    MPI_Allreduce(&local_max, &global_value, 1, MPI_DOUBLE, MPI_MAX, Chi::mpi.comm);
    op_value_ = global_value;
  }
}