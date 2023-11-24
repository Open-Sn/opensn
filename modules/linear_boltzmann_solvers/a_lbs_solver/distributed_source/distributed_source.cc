#include "distributed_source.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/functions/spatial_material_function.h"

#include "framework/runtime.h"
#include "framework/object_factory.h"

namespace opensn::lbs
{

OpenSnRegisterObject(lbs, DistributedSource);

InputParameters
DistributedSource::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("General implementation of an arbitrary multi-group "
                               "distributed source. A distributed source is given by a "
                               "logical volume and a Lua function handle. Currently, only "
                               "isotropic distributed sources are allowed.");
  params.SetClassName("Distributed Source");
  params.SetDocGroup("LBSUtilities");

  params.AddRequiredParameter<size_t>("logical_volume_handle",
                                      "Handle to the logical volume the distributed source is "
                                      "defined within.");
  params.AddOptionalParameter(
    "function_handle",
    SIZE_T_INVALID,
    "Handle to a ResponseFunction object to be used to define the source.");

  return params;
}

DistributedSource::DistributedSource(const InputParameters& params)
  : logical_volume_ptr_(GetStackItemPtrAsType<opensn::LogicalVolume>(
      object_stack, params.GetParamValue<size_t>("logical_volume_handle"))),
    function_(params.ParametersAtAssignment().Has("function_handle")
                ? GetStackItemPtrAsType<opensn::SpatialMaterialFunction>(
                    object_stack, params.GetParamValue<size_t>("function_handle"))
                : nullptr)
{
}

void
lbs::DistributedSource::Initialize(const LBSSolver& lbs_solver)
{
  subscribers_.clear();
  for (const auto& cell : lbs_solver.Grid().local_cells)
    if (logical_volume_ptr_->Inside(cell.centroid_)) subscribers_.push_back(cell.local_id_);

  num_local_subsribers_ = subscribers_.size();
  mpi_comm.all_reduce(num_local_subsribers_, num_global_subscribers_, mpi::op::sum<size_t>());

  log.LogAll() << "Distributed source has " << num_local_subsribers_
               << " subscribing cells on processor " << opensn::mpi_comm.rank() << ".";
  log.Log() << "Distributed source has " << num_global_subscribers_ << " total subscribing cells.";
}

std::vector<double>
lbs::DistributedSource::operator()(const Cell& cell, const Vector3& xyz, const int num_groups) const
{
  if (std::count(subscribers_.begin(), subscribers_.end(), cell.local_id_) == 0)
    return std::vector<double>(num_groups, 0.0);
  else if (not function_)
    return std::vector<double>(num_groups, 1.0);
  else
    return function_->Evaluate(xyz, cell.material_id_, num_groups);
}

} // namespace opensn::lbs
