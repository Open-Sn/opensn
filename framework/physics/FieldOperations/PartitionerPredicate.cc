#include "PartitionerPredicate.h"

#include "graphs/GraphPartitioner.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiObjectFactory.h"


namespace chi_physics::field_operations
{

RegisterChiObject(chi_physics::field_operations, PartitionerPredicate);

chi::InputParameters PartitionerPredicate::GetInputParameters()
{
  chi::InputParameters params = FieldOperation::GetInputParameters();

  params.SetGeneralDescription("Field operation that will write to the "
                               "field the result of a partitioning operation");
  params.SetDocGroup("DocFieldOperation");

  params.AddRequiredParameter<size_t>(
    "partitioner",
    "Handle to a GraphPartitioner object to use for parallel partitioning.");

  params.AddRequiredParameter<size_t>(
    "result_field",
    "Handle to, or name of, the field function that should "
    "receive the result of the operation.");
  params.SetParameterTypeMismatchAllowed("result_field");

  params.AddRequiredParameter<size_t>(
    "num_partitions", "Number of parts to apply to the partitioning");

  params.AddOptionalParameter(
    "result_component",
    0,
    "Resulting component into which the result will be written.");

  return params;
}

PartitionerPredicate::PartitionerPredicate(const chi::InputParameters& params)
  : FieldOperation(params),
    partitioner_(Chi::GetStackItem<chi::GraphPartitioner>(
      Chi::object_stack,
      params.GetParamValue<size_t>("partitioner"),
      __FUNCTION__)),
    result_field_param_(params.GetParam("result_field")),
    num_partitions_(params.GetParamValue<size_t>("num_partitions")),
    result_component_(params.GetParamValue<size_t>("result_component"))
{
}

void PartitionerPredicate::Execute()
{
  std::shared_ptr<FieldFunctionGridBased> grid_ff_ptr;
  if (result_field_param_.Type() == chi::ParameterBlockType::INTEGER)
  {
    const size_t handle = result_field_param_.GetValue<size_t>();
    auto ff_ptr = Chi::GetStackItemPtrAsType<FieldFunction>(
      Chi::field_function_stack, handle, __FUNCTION__);
    grid_ff_ptr = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_ptr);

    ChiLogicalErrorIf(
      not grid_ff_ptr,
      "Could not cast field function to FieldFunctionGridBased");
  }
  else if (result_field_param_.Type() == chi::ParameterBlockType::STRING)
  {
    const std::string ff_name = result_field_param_.GetValue<std::string>();
    for (auto& ff_ptr : Chi::field_function_stack)
      if (ff_ptr->TextName() == ff_name)
      {
        grid_ff_ptr = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_ptr);

        ChiLogicalErrorIf(
          not grid_ff_ptr,
          "Could not cast field function to FieldFunctionGridBased");
      }
  }

  ChiInvalidArgumentIf(
    not grid_ff_ptr, "Could not find the associated resulting field function");

  //============================================= Build cell graph and centroids
  const auto& sdm = grid_ff_ptr->GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  typedef std::vector<uint64_t> CellGraphNode;
  typedef std::vector<CellGraphNode> CellGraph;
  CellGraph cell_graph;
  std::vector<chi_mesh::Vector3> cell_centroids;

  const size_t num_local_cells = grid.local_cells.size();

  cell_graph.reserve(num_local_cells);
  cell_centroids.reserve(num_local_cells);
  {
    for (const auto& cell : grid.local_cells)
    {
      CellGraphNode cell_graph_node; // <-- Note A
      for (const auto& face : cell.faces_)
        if (face.has_neighbor_) cell_graph_node.push_back(face.neighbor_id_);

      cell_graph.push_back(cell_graph_node);
      cell_centroids.push_back(cell.centroid_);
    }
  }

  //============================================= Create partition
  auto cell_pids = partitioner_.Partition(
    cell_graph, cell_centroids, static_cast<int>(num_partitions_));

  auto& local_data = grid_ff_ptr->FieldVector();
  const auto& uk_man = grid_ff_ptr->GetUnknownManager();

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const int64_t cell_pid = cell_pids.at(cell.local_id_);

    for (uint32_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dof_map =
        sdm.MapDOFLocal(cell, i, uk_man, 0, result_component_);

      local_data[dof_map] = static_cast<double>(cell_pid);
    }
  } // for cell
}

} // namespace chi_physics::field_operations