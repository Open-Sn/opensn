// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/vacuum_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/isotropic_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/arbitrary_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/tgdsa.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cmath>
#include <list>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

namespace opensn
{

namespace
{

std::set<std::uint64_t>
GetGlobalUniqueBoundaryIDs(const std::shared_ptr<MeshContinuum>& grid, mpi::Communicator& mpi_comm)
{
  std::set<std::uint64_t> local_unique_bids_set;
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& face : cell->faces)
      if (not face.has_neighbor)
        local_unique_bids_set.insert(face.neighbor_id);

  const std::vector<std::uint64_t> local_unique_bids(local_unique_bids_set.begin(),
                                                     local_unique_bids_set.end());

  std::vector<std::uint64_t> recvbuf;
  mpi_comm.all_gather(local_unique_bids, recvbuf);

  std::set<std::uint64_t> global_unique_bids_set = local_unique_bids_set;
  global_unique_bids_set.insert(recvbuf.begin(), recvbuf.end());

  return global_unique_bids_set;
}

// RecursiveAngleSort's traversal order determines, for each reflected pair of angle sets, which
// set's angle indices get propagated to the other. With a single reflecting boundary this is
// order-independent but with opposing reflecting boundaries, the final angle-to-angleset
// assignment depends on visitation order. unsorted/sorted MUST order by a value stable across
// MPI ranks and repeated runs. For this purpose, we use AngleSet::GetID().
struct AngleSetIdLess
{
  bool operator()(const AngleSet* a, const AngleSet* b) const { return a->GetID() < b->GetID(); }
};

void
RecursiveAngleSort(AngleSet* angleset,
                   AngleAggregation& angle_agg,
                   const std::vector<std::vector<std::uint32_t>>& reflected_maps,
                   std::set<AngleSet*, AngleSetIdLess>& unsorted,
                   std::set<AngleSet*, AngleSetIdLess>& sorted)
{
  sorted.insert(angleset);
  std::uint32_t angle_zero = angleset->GetAngleIndices()[0];
  for (const auto& reflected_map : reflected_maps)
  {
    auto reflected_angle_zero = reflected_map[angle_zero];
    auto* reflected_angleset = angle_agg.GetAngleSetForAngleIndex(reflected_angle_zero);
    if (sorted.contains(reflected_angleset))
    {
      bool is_coherent = std::equal(angleset->GetAngleIndices().begin(),
                                    angleset->GetAngleIndices().end(),
                                    reflected_angleset->GetAngleIndices().begin(),
                                    [&](std::uint32_t angle, std::uint32_t reflected)
                                    { return reflected_map[angle] == reflected; });
      if (not is_coherent)
        throw std::logic_error("Cannot find an unanimous sort order for the angle set.\n");
    }
    else
    {
      std::transform(angleset->GetAngleIndices().begin(),
                     angleset->GetAngleIndices().end(),
                     reflected_angleset->GetAngleIndices().begin(),
                     [&](std::uint32_t angle) { return reflected_map[angle]; });
      reflected_angleset->SyncDeviceAngleIndices();
      unsorted.erase(reflected_angleset);
      RecursiveAngleSort(reflected_angleset, angle_agg, reflected_maps, unsorted, sorted);
    }
  }
}

} // namespace

const std::map<uint64_t, std::shared_ptr<SweepBoundary>>&
DiscreteOrdinatesProblem::GetSweepBoundaries() const
{
  return sweep_boundaries_;
}

const std::map<uint64_t, BoundaryDefinition>&
DiscreteOrdinatesProblem::GetBoundaryDefinitions() const
{
  return boundary_definitions_;
}

void
DiscreteOrdinatesProblem::SetBoundaryOptions(const std::vector<InputParameters>& boundary_params,
                                             bool clear_existing)
{
  if (clear_existing)
    boundary_definitions_.clear();

  for (const auto& params : boundary_params)
    UpdateBoundaryDefinition(params);

  if (clear_existing or not boundary_params.empty())
  {
    InitializeBoundaries();
    RebuildBoundaryRuntimeData();
  }
}

void
DiscreteOrdinatesProblem::UpdateBoundaryDefinition(const InputParameters& params)
{
  const auto boundary_name = params.GetParamValue<std::string>("name");
  const auto coord_sys = grid_->GetCoordinateSystem();
  const auto bnd_name_map = grid_->GetBoundaryNameMap();
  const auto mesh_type = grid_->GetType();

  // If we're using RZ, the user should use rmin/rmax/zmin/zmax and we'll
  // map internally to xmin/xmax/ymin/max
  std::string lookup_name = boundary_name;
  if (coord_sys == CoordinateSystemType::CYLINDRICAL and mesh_type == MeshType::ORTHOGONAL and
      grid_->GetDimension() == 2)
  {
    if (boundary_name != "rmin" and boundary_name != "rmax" and boundary_name != "zmin" and
        boundary_name != "zmax")
    {
      throw std::runtime_error(GetName() + ": Boundary name '" + boundary_name +
                               "' is invalid for cylindrical orthogonal meshes. "
                               "Use rmin, rmax, zmin, zmax.");
    }

    const std::map<std::string, std::string> rz_map = {
      {"rmin", "xmin"}, {"rmax", "xmax"}, {"zmin", "ymin"}, {"zmax", "ymax"}};
    const auto rz_it = rz_map.find(boundary_name);
    if (rz_it != rz_map.end())
      lookup_name = rz_it->second;
  }

  const auto it = bnd_name_map.find(lookup_name);
  if (it == bnd_name_map.end())
  {
    throw std::runtime_error("Boundary name \"" + boundary_name + "\" not found in mesh.");
  }
  const auto bid = it->second;
  boundary_definitions_[bid] = BoundaryDefinition(params, GetNumGroups());
}

void
DiscreteOrdinatesProblem::ClearBoundaries()
{
  boundary_definitions_.clear();
  InitializeBoundaries();
  RebuildBoundaryRuntimeData();
}

void
DiscreteOrdinatesProblem::RebuildBoundaryRuntimeData()
{
  if (boundary_runtime_data_initialized_)
    return;

  const bool solver_schemes_initialized =
    ags_solver_ or (not wgs_solvers_.empty()) or (not wgs_contexts_.empty());

  for (auto& [bid, boundary] : sweep_boundaries_)
    boundary->InitializeReflectingMap(groupsets_);
  if (use_gpus_)
    SortAngleSetsAngleIndices();
  for (auto& [bid, boundary] : sweep_boundaries_)
    boundary->InitializeAngleDependent(groupsets_);
  boundary_bank_.ShrinkToFit();
  boundary_bank_.DisableAllocation();
  InitializeBoundaryCarrier();
  if (sweep_type_ == "AAH" and use_gpus_)
    UpdateAAHD_FLUDSCommonDataWithBoundary();
  for (auto& groupset : groupsets_)
    groupset.angle_agg->SetupAngleSetDependencies();

  if (solver_schemes_initialized)
  {
    for (auto& groupset : groupsets_)
    {
      WGDSA::CleanUp(groupset);
      TGDSA::CleanUp(groupset);
      WGDSA::Init(*this, groupset);
      TGDSA::Init(*this, groupset);
    }
    ReinitializeSolverSchemes();
  }

  boundary_runtime_data_initialized_ = true;
}

Vector3
DiscreteOrdinatesProblem::ComputeReflectingBoundaryNormal(uint64_t bid) const
{
  const double EPSILON = 1.0e-12;
  std::unique_ptr<Vector3> n_ptr = nullptr;
  for (const auto& cell : grid_->GetLocalCells())
  {
    for (const auto& face : cell->faces)
    {
      if (not face.has_neighbor and face.neighbor_id == bid)
      {
        if (not n_ptr)
          n_ptr = std::make_unique<Vector3>(face.normal);
        if (std::fabs(face.normal.Dot(*n_ptr) - 1.0) > EPSILON)
          throw std::logic_error(GetName() +
                                 ": Not all face normals are, within tolerance, locally the same "
                                 "for the reflecting boundary condition requested");
      }
    }
  }

  const int local_has_bid = n_ptr != nullptr ? 1 : 0;
  const Vector3 local_normal = local_has_bid ? *n_ptr : Vector3(0.0, 0.0, 0.0);

  std::vector<int> locJ_has_bid(opensn::mpi_comm.size(), 1);
  std::vector<double> locJ_n_val(opensn::mpi_comm.size() * 3L, 0.0);

  mpi_comm.all_gather(local_has_bid, locJ_has_bid);
  std::vector<double> lnv = {local_normal.x, local_normal.y, local_normal.z};
  mpi_comm.all_gather(lnv.data(), 3, locJ_n_val.data(), 3);

  Vector3 global_normal;
  for (int j = 0; j < opensn::mpi_comm.size(); ++j)
  {
    if (locJ_has_bid[j])
    {
      int offset = 3 * j;
      const double* n = &locJ_n_val[offset];
      const Vector3 locJ_normal(n[0], n[1], n[2]);

      if (local_has_bid)
        if (std::fabs(local_normal.Dot(locJ_normal) - 1.0) > EPSILON)
          throw std::logic_error(GetName() +
                                 ": Not all face normals are, within tolerance, globally the same "
                                 "for the reflecting boundary condition requested");

      global_normal = locJ_normal;
    }
  }

  return global_normal;
}

void
DiscreteOrdinatesProblem::InitializeBoundaries()
{
  CALI_CXX_MARK_SCOPE("Boundaries");

  ResetBoundaryCarrier();
  boundary_runtime_data_initialized_ = false;
  has_reflecting_boundaries_ = false;
  has_time_dependent_boundaries_ = false;

  ValidateBoundaryConfiguration();

  // Determine boundary-ids involved in the problem
  const auto unique_bids_set = GetGlobalUniqueBoundaryIDs(grid_, mpi_comm);
  for (const auto bid : unique_bids_set)
  {
    if (boundary_definitions_.find(bid) == boundary_definitions_.end())
      boundary_definitions_.emplace(bid, BoundaryDefinition());
  }
  boundary_bank_ = BoundaryBank(groupsets_);

  sweep_boundaries_.clear();
  const auto coord_sys = MapGeometryTypeToCoordSys(geometry_type_);
  for (auto bid : unique_bids_set)
  {
    auto& bndry_def = boundary_definitions_.find(bid)->second;

    switch (bndry_def.type)
    {
      case LBSBoundaryType::VACUUM:
      {
        sweep_boundaries_[bid] = std::make_shared<VacuumBoundary>(boundary_bank_);
        break;
      }
      case LBSBoundaryType::ISOTROPIC:
      {
        sweep_boundaries_[bid] = std::make_shared<IsotropicBoundary>(boundary_bank_,
                                                                     groupsets_,
                                                                     bndry_def.group_strength,
                                                                     bndry_def.start_time,
                                                                     bndry_def.end_time);
        has_time_dependent_boundaries_ = true;
        break;
      }
      case LBSBoundaryType::ARBITRARY:
      {
        sweep_boundaries_[bid] = std::make_shared<ArbitraryBoundary>(
          boundary_bank_, groupsets_, bndry_def.time_angular_flux_function);
        has_time_dependent_boundaries_ = true;
        break;
      }
      case LBSBoundaryType::REFLECTING:
      {
        const auto global_normal = ComputeReflectingBoundaryNormal(bid);
        sweep_boundaries_[bid] = std::make_shared<ReflectingBoundary>(
          boundary_bank_, bid, grid_, groupsets_, global_normal, coord_sys);
        has_reflecting_boundaries_ = true;
        break;
      }
      default:
      {
        throw std::logic_error("Boundary type not implemented.");
      }
    }
  }

  for (auto& [bid, bndry] : sweep_boundaries_)
  {
    bndry->SetOpposingReflected(bid, sweep_boundaries_);
  }
}

void
DiscreteOrdinatesProblem::SortAngleSetsAngleIndices()
{
  for (auto& groupset : groupsets_)
  {
    // build the total reflected map
    std::vector<std::vector<std::uint32_t>> reflected_maps;
    for (auto& [bid, boundary] : sweep_boundaries_)
      boundary->GetReflectedMap(groupset.id, reflected_maps);

    // check for each map and angle set that the mapped angle indices is also another angleset
    auto& angle_agg = *(groupset.angle_agg);
    std::list<std::set<std::uint32_t>> angle_indices_list;
    for (const auto& angleset : angle_agg)
    {
      const auto& angle_indices = angleset->GetAngleIndices();
      angle_indices_list.emplace_back(angle_indices.begin(), angle_indices.end());
    }
    for (const auto& reflected_map : reflected_maps)
    {
      std::list<std::set<std::uint32_t>> angle_indices_list_copy(angle_indices_list);
      for (auto it_in = angle_indices_list_copy.begin(); it_in != angle_indices_list_copy.end();)
      {
        std::set<std::uint32_t> reflected;
        std::transform(it_in->begin(),
                       it_in->end(),
                       std::inserter(reflected, reflected.end()),
                       [&](std::uint32_t n) { return reflected_map[n]; });
        auto it_out =
          std::find(angle_indices_list_copy.begin(), angle_indices_list_copy.end(), reflected);
        if (it_out == it_in or it_out == angle_indices_list_copy.end())
          throw std::logic_error("Angleset parity was broken for one of the reflected boundaries.");
        angle_indices_list_copy.erase(it_out);
        it_in = angle_indices_list_copy.erase(it_in);
      }
    }
    angle_indices_list.clear();

    // sort angle indices in anglesets. Ordering is by AngleSetIdLess (GetID()) so this
    // recursion is deterministic and consistent across MPI ranks
    std::set<AngleSet*, AngleSetIdLess> sorted_anglesets, unsorted_anglesets;
    std::transform(angle_agg.begin(),
                   angle_agg.end(),
                   std::inserter(unsorted_anglesets, unsorted_anglesets.end()),
                   [](auto& angleset) { return angleset.get(); });
    while (not unsorted_anglesets.empty())
    {
      AngleSet* angleset = *unsorted_anglesets.begin();
      unsorted_anglesets.erase(unsorted_anglesets.begin());
      RecursiveAngleSort(angleset, angle_agg, reflected_maps, unsorted_anglesets, sorted_anglesets);
    }
  }
}

#ifndef __OPENSN_WITH_GPU__
void
DiscreteOrdinatesProblem::InitializeBoundaryCarrier()
{
}

void
DiscreteOrdinatesProblem::TransferDeviceBoundaryData(int groupset_id,
                                                     bool host_to_device,
                                                     bool force)
{
}

void
DiscreteOrdinatesProblem::ResetBoundaryCarrier()
{
}

void
DiscreteOrdinatesProblem::UpdateAAHD_FLUDSCommonDataWithBoundary()
{
  throw std::runtime_error("DiscreteOrdinatesProblem::UpdateAAHD_FLUDSCommonDataWithBoundary : "
                           "OPENSN_WITH_CUDA not enabled.");
}
#endif

} // namespace opensn
