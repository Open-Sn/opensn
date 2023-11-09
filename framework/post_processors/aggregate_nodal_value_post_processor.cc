#include "framework/post_processors/aggregate_nodal_value_post_processor.h"

#include "framework/object_factory.h"

#include "framework/physics/field_function/field_function_grid_based.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/event_system/event.h"

namespace chi
{

RegisterChiObject(chi, AggregateNodalValuePostProcessor);

InputParameters
AggregateNodalValuePostProcessor::GetInputParameters()
{
  InputParameters params = PostProcessor::GetInputParameters();
  params += chi_physics::GridBasedFieldFunctionInterface::GetInputParameters();
  params += chi_mesh::LogicalVolumeInterface::GetInputParameters();

  params.SetGeneralDescription("Gets the max/min/avg nodal value of a field function "
                               "among nodal values.");
  params.SetDocGroup("doc_PostProcessors");

  params.AddRequiredParameter<std::string>("operation", "The required operation to be performed.");

  using namespace chi_data_types;
  params.ConstrainParameterRange("operation", AllowableRangeList::New({"max", "min", "avg"}));

  return params;
}

AggregateNodalValuePostProcessor::AggregateNodalValuePostProcessor(const InputParameters& params)
  : PostProcessor(params, PPType::SCALAR),
    chi_physics::GridBasedFieldFunctionInterface(params),
    chi_mesh::LogicalVolumeInterface(params),
    operation_(params.GetParamValue<std::string>("operation"))
{
}

void
AggregateNodalValuePostProcessor::Initialize()
{
  const auto* grid_field_function = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not grid_field_function,
                    "Attempted to access invalid field"
                    "function");

  const auto& grid = grid_field_function->GetSpatialDiscretization().Grid();

  const auto* logical_volume_ptr_ = GetLogicalVolume();
  if (logical_volume_ptr_ == nullptr)
  {
    cell_local_ids_.reserve(grid.local_cells.size());
    for (const auto& cell : grid.local_cells)
      cell_local_ids_.push_back(cell.local_id_);
  }
  else
  {
    for (const auto& cell : grid.local_cells)
      if (logical_volume_ptr_->Inside(cell.centroid_)) cell_local_ids_.push_back(cell.local_id_);
  }

  initialized_ = true;
}

void
AggregateNodalValuePostProcessor::Execute(const Event& event_context)
{
  if (not initialized_) Initialize();

  const auto* grid_field_function = GetGridBasedFieldFunction();

  ChiLogicalErrorIf(not grid_field_function,
                    "Attempted to access invalid field"
                    "function");

  const auto& ref_ff = *grid_field_function;
  const auto& sdm = ref_ff.GetSpatialDiscretization();
  const auto& grid = sdm.Grid();

  const auto& uk_man = ref_ff.GetUnknownManager();
  const auto uid = 0;
  const auto cid = 0;

  const auto field_data = ref_ff.GetGhostedFieldVector();
  const size_t num_local_dofs =
    ref_ff.GetSpatialDiscretization().GetNumLocalDOFs(ref_ff.GetUnknownManager());

  double local_max_value = 0.0;
  double local_min_value = 0.0;
  double local_accumulation = 0.0;
  bool first_local = true;
  for (const uint64_t cell_local_id : cell_local_ids_)
  {
    const auto& cell = grid.local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
      if (imap >= 0 and imap < num_local_dofs)
      {
        const double field_value = field_data[imap];
        if (first_local)
        {
          local_max_value = field_value;
          local_min_value = field_value;
          first_local = false;
        }

        local_max_value = std::max(local_max_value, field_value);
        local_min_value = std::min(local_min_value, field_value);
        local_accumulation += field_value;
      }
    } // for i
  }   // for cell-id

  if (operation_ == "max")
  {
    double globl_max_value;
    MPI_Allreduce(&local_max_value, &globl_max_value, 1, MPI_DOUBLE, MPI_MAX, Chi::mpi.comm);

    value_ = ParameterBlock("", globl_max_value);
  }
  else if (operation_ == "min")
  {
    double globl_min_value;
    MPI_Allreduce(&local_min_value, &globl_min_value, 1, MPI_DOUBLE, MPI_MIN, Chi::mpi.comm);

    value_ = ParameterBlock("", globl_min_value);
  }
  else if (operation_ == "avg")
  {
    double globl_accumulation;
    MPI_Allreduce(&local_accumulation, &globl_accumulation, 1, MPI_DOUBLE, MPI_SUM, Chi::mpi.comm);

    const size_t num_globl_dofs =
      ref_ff.GetSpatialDiscretization().GetNumGlobalDOFs(ref_ff.GetUnknownManager());
    value_ = ParameterBlock("", globl_accumulation / double(num_globl_dofs));
  }
  else
    ChiLogicalError("Unsupported operation type \"" + operation_ + "\".");

  const int event_code = event_context.Code();
  if (event_code == Event::SolverInitialized or event_code == Event::SolverAdvanced)
  {
    const auto& event_params = event_context.Parameters();

    if (event_params.Has("timestep_index") and event_params.Has("time"))
    {
      const size_t index = event_params.GetParamValue<size_t>("timestep_index");
      const double time = event_params.GetParamValue<double>("time");
      TimeHistoryEntry entry{index, time, value_};
      time_history_.push_back(std::move(entry));
    }
  }
}

} // namespace chi
