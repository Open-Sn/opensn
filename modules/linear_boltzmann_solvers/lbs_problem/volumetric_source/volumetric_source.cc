// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/volumetric_source/volumetric_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/math/functions/function.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/runtime.h"
#include "framework/object_factory.h"
#include <memory>
#include <limits>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, VolumetricSource);

int VolumetricSource::next_id_ = 0;

InputParameters
VolumetricSource::GetInputParameters()
{
  InputParameters params;

  params.SetGeneralDescription("General implementation of an arbitrary multi-group "
                               "volumetric source. Currently, only isotropic volumetric "
                               "sources are allowed.");
  params.SetClassName("Volumetric Source");

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
  params.AddOptionalParameter<std::shared_ptr<VectorSpatialFunction>>(
    "func",
    std::shared_ptr<VectorSpatialFunction>{},
    "SpatialMaterialFunction object to be used to define the source.");
  params.AddOptionalParameter("start_time",
                              -std::numeric_limits<double>::infinity(),
                              "Time at which the source becomes active.");
  params.AddOptionalParameter(
    "end_time", std::numeric_limits<double>::infinity(), "Time at which the source is inactive.");

  return params;
}

std::shared_ptr<VolumetricSource>
VolumetricSource::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<VolumetricSource>("lbs::VolumetricSource", params);
}

VolumetricSource::VolumetricSource(const InputParameters& params)
  : id_(next_id_++),
    block_ids_(params.GetParamVectorValue<unsigned int>("block_ids")),
    logvol_(params.GetSharedPtrParam<LogicalVolume>("logical_volume", false)),
    strength_(params.GetParamVectorValue<double>("group_strength")),
    function_(params.GetSharedPtrParam<VectorSpatialFunction>("func", false)),
    start_time_(params.GetParamValue<double>("start_time")),
    end_time_(params.GetParamValue<double>("end_time"))
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
VolumetricSource::Initialize(const LBSProblem& lbs_problem)
{
  // Set the source strength vector
  if (not function_ and not strength_.empty())
    if (strength_.size() != lbs_problem.GetNumGroups())
      throw std::invalid_argument("The number of groups in the source strength vector must "
                                  "match the number of groups in the solver the source is "
                                  "attached to.");

  // Set cell subscribers based on logical volumes, block IDs, or both
  subscribers_.clear();
  if (logvol_ and block_ids_.empty())
  {
    std::set<unsigned int> blk_ids;
    for (const auto& cell : lbs_problem.GetGrid()->local_cells)
      if (logvol_->Inside(cell.centroid))
      {
        blk_ids.insert(cell.block_id);
        subscribers_.push_back(cell.local_id);
      }
  }
  else if (not logvol_ and not block_ids_.empty())
  {
    for (const auto& cell : lbs_problem.GetGrid()->local_cells)
      if (std::find(block_ids_.begin(), block_ids_.end(), cell.block_id) != block_ids_.end())
        subscribers_.push_back(cell.local_id);
  }
  else
  {
    for (const auto& cell : lbs_problem.GetGrid()->local_cells)
      if (logvol_->Inside(cell.centroid) and
          std::find(block_ids_.begin(), block_ids_.end(), cell.block_id) != block_ids_.end())
        subscribers_.push_back(cell.local_id);
  }

  num_local_subsribers_ = subscribers_.size();
  mpi_comm.all_reduce(num_local_subsribers_, num_global_subscribers_, mpi::op::sum<size_t>());

  log.LogAllVerbose1() << "Volumetric source #" << id_ << " has " << num_local_subsribers_
                       << " subscribing cells on processor " << opensn::mpi_comm.rank() << ".";
  log.Log() << "Volumetric source #" << id_ << " has " << num_global_subscribers_
            << " total subscribing cells.";
}

std::vector<double>
VolumetricSource::operator()(const Cell& cell,
                             const Vector3& xyz,
                             const std::size_t num_groups) const
{
  if (std::count(subscribers_.begin(), subscribers_.end(), cell.local_id) == 0)
    return std::vector<double>(num_groups, 0.0); // NOLINT
  else if (not function_)
    return strength_;
  else
    return (*function_)(xyz, num_groups);
}

bool
VolumetricSource::IsActive(double time) const
{
  return time >= start_time_ && time <= end_time_;
}

} // namespace opensn
