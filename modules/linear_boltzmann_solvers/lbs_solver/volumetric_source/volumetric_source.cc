// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/volumetric_source/volumetric_source.h"
#include "framework/math/functions/function.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/functions/vector_spatial_function.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"
#include <memory>

namespace opensn
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
  params.AddOptionalParameterArray("group_strength",
                                   std::vector<double>(),
                                   "An array of multi-group source strength values. Note that this "
                                   "is only used when a function is not provided.");
  params.AddOptionalParameter<std::shared_ptr<LogicalVolume>>(
    "logical_volume",
    std::shared_ptr<LogicalVolume>{},
    "Handle to the logical volume the volumetric source is defined within.");
  params.AddOptionalParameter<std::shared_ptr<Function>>(
    "func",
    std::shared_ptr<Function>{},
    "SpatialMaterialFunction object to be used to define the source.");

  return params;
}

std::shared_ptr<VolumetricSource>
VolumetricSource::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<VolumetricSource>("lbs::VolumetricSource", params);
}

VolumetricSource::VolumetricSource(const InputParameters& params)
  : block_ids_(params.GetParamVectorValue<int>("block_ids")),
    logvol_(params.GetParamValue<std::shared_ptr<LogicalVolume>>("logical_volume")),
    strength_(params.GetParamVectorValue<double>("group_strength")),
    function_(std::dynamic_pointer_cast<VectorSpatialFunction>(
      params.GetParamValue<std::shared_ptr<Function>>("func")))
{
  if (not logvol_ and block_ids_.empty())
    throw std::invalid_argument("A volumetric source must be defined with a logical volume, "
                                "block IDs, or both. Neither were specified.");
  if (function_ and not strength_.empty())
    throw std::invalid_argument(
      "If a function is provided, the source strength should not be set.");
  if (not function_ and strength_.empty())
    throw std::invalid_argument("Either a function or the source strength must be provided.");
}

void
VolumetricSource::Initialize(const LBSSolver& lbs_solver)
{
  // Set the source strength vector
  if (not function_ and not strength_.empty())
    if (strength_.size() != lbs_solver.GetNumGroups())
      throw std::invalid_argument("The number of groups in the source strength vector must "
                                  "match the number of groups in the solver the source is "
                                  "attached to.");

  // Set cell subscribers based on logical volumes, block IDs, or both
  subscribers_.clear();
  if (logvol_ and block_ids_.empty())
  {
    std::set<int> blk_ids;
    for (const auto& cell : lbs_solver.GetGrid().local_cells)
      if (logvol_->Inside(cell.centroid))
      {
        blk_ids.insert(cell.material_id);
        subscribers_.push_back(cell.local_id);
      }
  }
  else if (not logvol_ and not block_ids_.empty())
  {
    for (const auto& cell : lbs_solver.GetGrid().local_cells)
      if (std::find(block_ids_.begin(), block_ids_.end(), cell.material_id) != block_ids_.end())
        subscribers_.push_back(cell.local_id);
  }
  else
  {
    for (const auto& cell : lbs_solver.GetGrid().local_cells)
      if (logvol_->Inside(cell.centroid) and
          std::find(block_ids_.begin(), block_ids_.end(), cell.material_id) != block_ids_.end())
        subscribers_.push_back(cell.local_id);
  }

  num_local_subsribers_ = subscribers_.size();
  mpi_comm.all_reduce(num_local_subsribers_, num_global_subscribers_, mpi::op::sum<size_t>());

  log.LogAll() << "Volumetric source has " << num_local_subsribers_
               << " subscribing cells on processor " << opensn::mpi_comm.rank() << ".";
  log.Log() << "Volumetric source has " << num_global_subscribers_ << " total subscribing cells.";
}

std::vector<double>
VolumetricSource::operator()(const Cell& cell, const Vector3& xyz, const int num_groups) const
{
  if (std::count(subscribers_.begin(), subscribers_.end(), cell.local_id) == 0)
    return std::vector<double>(num_groups, 0.0);
  else if (not function_)
    return strength_;
  else
    return function_->Evaluate(xyz, num_groups);
}

} // namespace opensn
