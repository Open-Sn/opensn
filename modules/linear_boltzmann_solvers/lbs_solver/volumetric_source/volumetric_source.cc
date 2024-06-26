// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/functions/spatial_material_function.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"

namespace opensn::lbs
{

OpenSnRegisterObjectInNamespace(lbs, VolumetricSource);

InputParameters
VolumetricSource::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("General implementation of an arbitrary multi-group "
                               "volumetric source. A volumetric source is given by a "
                               "logical volume and a Lua function handle. Currently, only "
                               "isotropic volumetric sources are allowed.");
  params.SetClassName("Volumetric Source");
  params.SetDocGroup("LBSUtilities");

  params.AddOptionalParameterArray(
    "block_ids",
    std::vector<int>(),
    "An array of block IDs the volumetric source is present within.");
  params.AddOptionalParameter(
    "logical_volume_handle",
    SIZE_T_INVALID,
    "Handle to the logical volume the volumetric source is defined within.");
  params.AddOptionalParameter(
    "function_handle",
    SIZE_T_INVALID,
    "Handle to a SpatialMaterialFunction object to be used to define the source.");

  return params;
}

VolumetricSource::VolumetricSource(const InputParameters& params)
  : block_ids_(params.GetParamVectorValue<int>("block_ids")),
    logvol_(params.ParametersAtAssignment().Has("logical_volume_handle")
              ? GetStackItemPtrAsType<LogicalVolume>(
                  object_stack, params.GetParamValue<size_t>("logical_volume_handle"))
              : nullptr),
    function_(params.ParametersAtAssignment().Has("function_handle")
                ? GetStackItemPtrAsType<SpatialMaterialFunction>(
                    object_stack, params.GetParamValue<size_t>("function_handle"))
                : nullptr)
{
  if (not logvol_ and block_ids_.empty())
    throw std::invalid_argument("A volumetric source must be defined with a logical volume, "
                                "block IDs, or both. Neither were specified.");
}

void
lbs::VolumetricSource::Initialize(const LBSSolver& lbs_solver)
{
  subscribers_.clear();

  // Subscribers from logical volume only
  if (logvol_ and block_ids_.empty())
  {
    std::set<int> blk_ids;
    for (const auto& cell : lbs_solver.Grid().local_cells)
      if (logvol_->Inside(cell.centroid_))
      {
        blk_ids.insert(cell.material_id_);
        subscribers_.push_back(cell.local_id_);
      }
  }

  // Subscribers from logical volume and block IDs
  else if (logvol_ and not block_ids_.empty())
  {
    for (const auto& cell : lbs_solver.Grid().local_cells)
      if (logvol_->Inside(cell.centroid_) and
          std::find(block_ids_.begin(), block_ids_.end(), cell.material_id_) != block_ids_.end())
        subscribers_.push_back(cell.local_id_);
  }

  // Subscribers from block IDs only.
  else
  {
    for (const auto& cell : lbs_solver.Grid().local_cells)
      if (std::find(block_ids_.begin(), block_ids_.end(), cell.material_id_) != block_ids_.end())
        subscribers_.push_back(cell.local_id_);
  }

  num_local_subsribers_ = subscribers_.size();
  mpi_comm.all_reduce(num_local_subsribers_, num_global_subscribers_, mpi::op::sum<size_t>());

  log.LogAll() << "Distributed source has " << num_local_subsribers_
               << " subscribing cells on processor " << opensn::mpi_comm.rank() << ".";
  log.Log() << "Distributed source has " << num_global_subscribers_ << " total subscribing cells.";
}

std::vector<double>
lbs::VolumetricSource::operator()(const Cell& cell, const Vector3& xyz, const int num_groups) const
{
  std::vector<double> src;
  if (std::count(subscribers_.begin(), subscribers_.end(), cell.local_id_) == 0)
    src.assign(num_groups, 0.0);
  else if (not function_)
    src.assign(num_groups, 1.0);
  else
    src = function_->Evaluate(xyz, cell.material_id_, num_groups);
  return src;
}

} // namespace opensn::lbs
