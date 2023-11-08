#include "framework/mesh/FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"
#include "framework/physics/FieldFunction/fieldfunction_gridbased.h"
#include "framework/math/SpatialDiscretization/SpatialDiscretization.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/math/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "framework/math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/console/chi_console.h"

namespace chi_mesh
{

void
FieldFunctionInterpolationVolume::Initialize()
{
  Chi::log.Log0Verbose1() << "Initializing volume interpolator.";
  // Check grid available
  if (field_functions_.empty())
    throw std::logic_error("Unassigned field function in volume field "
                           "function interpolator.");

  if (logical_volume_ == nullptr)
    throw std::logic_error("Unassigned logical volume in volume field function"
                           "interpolator.");

  const auto& grid = field_functions_.front()->GetSpatialDiscretization().Grid();

  // Find cells inside volume
  for (const auto& cell : grid.local_cells)
    if (logical_volume_->Inside(cell.centroid_))
      cell_local_ids_inside_logvol_.push_back(cell.local_id_);
}

void
FieldFunctionInterpolationVolume::Execute()
{
  const auto& ref_ff = *field_functions_.front();
  const auto& sdm = ref_ff.GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.GetUnknownManager();
  const auto uid = 0;
  const auto cid = ref_component_;

  using namespace ff_interpolation;
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
#ifdef OPENSN_WITH_LUA
      if (op_type_ >= Operation::OP_SUM_LUA and op_type_ <= Operation::OP_MAX_LUA)
        function_value = CallLuaFunction(ff_value, cell.material_id_);
#endif

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

#ifdef OPENSN_WITH_LUA
double
FieldFunctionInterpolationVolume::CallLuaFunction(double ff_value, int mat_id) const
{
  lua_State* L = Chi::console.GetConsoleState();
  double ret_val = 0.0;

  lua_getglobal(L, op_lua_func_.c_str());
  lua_pushnumber(L, ff_value);
  lua_pushnumber(L, mat_id);

  // 2 arguments, 1 result, 0=original error object
  if (lua_pcall(L, 2, 1, 0) == 0) { ret_val = lua_tonumber(L, -1); }
  lua_pop(L, 1);

  return ret_val;
}
#endif

} // namespace chi_mesh
