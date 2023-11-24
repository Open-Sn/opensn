#include "framework/mesh/field_function_interpolation/ffinter_volume.h"
#include "framework/physics/field_function/field_function_grid_based.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/vector_ghost_communicator/vector_ghost_communicator.h"
#include "framework/math/spatial_discretization/finite_element/quadrature_point_data.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

void
FieldFunctionInterpolationVolume::SetOperationFunction(
  std::shared_ptr<ScalarMaterialFunction> function)
{
  oper_function_ = function;
}

void
FieldFunctionInterpolationVolume::Initialize()
{
  log.Log0Verbose1() << "Initializing volume interpolator.";
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
      if (op_type_ >= FieldFunctionInterpolationOperation::OP_SUM_FUNC and
          op_type_ <= FieldFunctionInterpolationOperation::OP_MAX_FUNC)
        function_value = oper_function_->Evaluate(ff_value, cell.material_id_);

      local_volume += qp_data.JxW(qp);
      local_sum += function_value * qp_data.JxW(qp);
      local_max = std::fmax(ff_value, local_max);
      local_min = std::fmin(ff_value, local_min);
    } // for qp
  }   // for cell-id

  if (op_type_ == FieldFunctionInterpolationOperation::OP_SUM or
      op_type_ == FieldFunctionInterpolationOperation::OP_SUM_FUNC)
  {
    mpi_comm.all_reduce(local_sum, op_value_, mpi::op::sum<double>());
  }
  if (op_type_ == FieldFunctionInterpolationOperation::OP_AVG or
      op_type_ == FieldFunctionInterpolationOperation::OP_AVG_FUNC)
  {
    double local_data[] = {local_volume, local_sum};
    double global_data[] = {0.0, 0.0};

    mpi_comm.all_reduce(local_data, 2, global_data, mpi::op::sum<double>());
    double global_volume = global_data[0];
    double global_sum = global_data[1];
    op_value_ = global_sum / global_volume;
  }
  if (op_type_ == FieldFunctionInterpolationOperation::OP_MAX or
      op_type_ == FieldFunctionInterpolationOperation::OP_MAX_FUNC)
  {
    mpi_comm.all_reduce(local_max, op_value_, mpi::op::max<double>());
  }
}

} // namespace opensn
