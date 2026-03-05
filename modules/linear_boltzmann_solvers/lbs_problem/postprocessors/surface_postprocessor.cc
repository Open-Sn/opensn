// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/surface_postprocessor.h"
#include "framework/object_factory.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>

namespace opensn
{

namespace
{

auto
ComputeInflow(const LBSGroupset& groupset,
              const Cell& cell,
              const CellMapping& cell_mapping,
              SweepBoundary& bndry,
              const std::vector<Vector<double>>& IntS_shapeI,
              const CellFace& face,
              unsigned int face_idx,
              unsigned int g)
{
  double inflow = 0;
  for (int n = 0; n < groupset.quadrature->omegas.size(); ++n)
  {
    const auto& omega = groupset.quadrature->omegas[n];
    const double wt = groupset.quadrature->weights[n];
    const double mu = omega.Dot(face.normal);
    if (mu >= 0.)
      continue;

    for (int fi = 0; fi < face.vertex_ids.size(); ++fi)
    {
      const int i = cell_mapping.MapFaceNode(face_idx, fi);
      const auto& IntFi_shapeI = IntS_shapeI[face_idx](i);

      const double* psi_ptr = bndry.PsiIncoming(cell.local_id, face_idx, fi, n, g);
      const double psi = (psi_ptr != nullptr) ? *psi_ptr : 0.0;
      inflow += mu * wt * psi * IntFi_shapeI;
    }
  }
  return inflow;
};
} // namespace

OpenSnRegisterObjectInNamespace(lbs, SurfacePostprocessor);

InputParameters
SurfacePostprocessor::GetInputParameters()
{
  InputParameters params;
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "A handle to an existing LBS problem.");
  params.AddRequiredParameterArray(
    "boundaries",
    "Boundary restriction for the postprocessor. Empty/unspecified means no block restriction.");
  params.AddOptionalParameterArray("logical_volumes",
                                   std::vector<std::shared_ptr<LogicalVolume>>{},
                                   "Logical volume to restrict the computation to.");
  params.AddRequiredParameter<std::string>(
    "value_type", "Type of value to compute: 'integral', 'max', 'min', or 'avg'");
  params.AddRequiredParameter<std::string>(
    "current_type", "Type of value to compute: 'incoming', 'outgoing', 'net'");
  params.AddOptionalParameter(
    "group", 0, "Single group to compute (mutually exclusive with groupset).");
  params.AddOptionalParameter(
    "groupset", 0, "Single groupset to compute (mutually exclusive with group).");

  params.ConstrainParameterRange("value_type",
                                 AllowableRangeList::New({"integral", "max", "min", "avg"}));
  params.ConstrainParameterRange("current_type",
                                 AllowableRangeList::New({"incoming", "outgoing", "net"}));

  return params;
}

std::shared_ptr<SurfacePostprocessor>
SurfacePostprocessor::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SurfacePostprocessor>("lbs::SurfacePostprocessor", params);
}

SurfacePostprocessor::SurfacePostprocessor(const InputParameters& params)
  : do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem")),
    boundary_names_(params.GetParamVectorValue<std::string>("boundaries")),
    logical_volumes_(params.GetParamVectorValue<std::shared_ptr<LogicalVolume>>("logical_volumes")),
    selected_group_(params.IsParameterValid("group")
                      ? std::make_optional(params.GetParamValue<unsigned int>("group"))
                      : std::nullopt),
    selected_groupset_(params.IsParameterValid("groupset")
                         ? std::make_optional(params.GetParamValue<unsigned int>("groupset"))
                         : std::nullopt)
{
  const auto& grid = do_problem_->GetGrid();
  const auto& bnd_id_map = grid->GetBoundaryNameMap();
  for (auto& name : boundary_names_)
  {
    const auto& id = bnd_id_map.at(name);
    boundary_ids_.push_back(id);
  }
  // const auto & bnd_ids = do_problem_->GetSweepBoundaries();

  if (selected_group_.has_value() && selected_groupset_.has_value())
    throw std::invalid_argument("'group' and 'groupset' cannot be specified together");

  if (selected_group_.has_value() && selected_group_.value() >= do_problem_->GetNumGroups())
    throw std::invalid_argument("'group' must be less than " +
                                std::to_string(do_problem_->GetNumGroups()));

  if (selected_groupset_.has_value() &&
      selected_groupset_.value() >= do_problem_->GetGroupsets().size())
    throw std::invalid_argument("'groupset' must be less than " +
                                std::to_string(do_problem_->GetGroupsets().size()));

  const auto value_type_str = params.GetParamValue<std::string>("value_type");
  if (value_type_str == "max")
    value_type_ = ValueType::MAX;
  else if (value_type_str == "min")
    value_type_ = ValueType::MIN;
  else if (value_type_str == "integral")
    value_type_ = ValueType::INTEGRAL;
  else if (value_type_str == "avg")
    value_type_ = ValueType::AVERAGE;

  const auto current_type_str = params.GetParamValue<std::string>("current_type");
  if (current_type_str == "incoming")
    current_type_ = CurrentType::INCOMING;
  else if (current_type_str == "outgoing")
    current_type_ = CurrentType::OUTGOING;
  else if (current_type_str == "net")
    current_type_ = CurrentType::NET;

  CreateSpatialRestriction();
  CreateEnergyRestriction();

  auto n_lvs = std::max(static_cast<std::size_t>(1), logical_volumes_.size());
  values_.resize(n_lvs);
  for (auto& vals : values_)
    vals.resize(groups_.size());
}

void
SurfacePostprocessor::CreateSpatialRestriction()
{
  const auto& grid = do_problem_->GetGrid();
  const auto& groupsets = do_problem_->GetGroupsets();

  std::vector<SideFace> all_side_faces;
  for (const auto& cell : grid->local_cells)
  {
    for (int f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      // boundary face?
      if (not face.has_neighbor)
      {
        if (std::find(boundary_ids_.begin(), boundary_ids_.end(), face.neighbor_id) !=
            boundary_ids_.end())
          all_side_faces.push_back({cell.local_id, f});
      }
    }
  }

  if (logical_volumes_.empty())
  {
    side_faces_.resize(1);
    side_faces_[0] = all_side_faces;
  }
  else
  {
    side_faces_.resize(logical_volumes_.size());
    for (unsigned int i = 0; i < logical_volumes_.size(); ++i)
      side_faces_[i] = GetLogicalVolumeSides(all_side_faces, logical_volumes_[i]);
  }
}

std::vector<SurfacePostprocessor::SideFace>
SurfacePostprocessor::GetLogicalVolumeSides(const std::vector<SideFace>& all_side_faces,
                                            std::shared_ptr<LogicalVolume> log_vol)
{
  const auto& grid = do_problem_->GetGrid();
  std::vector<SideFace> side_faces;
  for (const auto& [cell, face_idx] : all_side_faces)
  {
    const auto& face = grid->local_cells[cell].faces[face_idx];
    if (log_vol->Inside(face.centroid))
      side_faces.emplace_back(cell, face_idx);
  }
  return side_faces;
}

void
SurfacePostprocessor::CreateEnergyRestriction()
{
  if (selected_group_.has_value())
  {
    const auto g = selected_group_.value();
    groups_.push_back(g);
    for (unsigned int gs_id = 0; gs_id < do_problem_->GetNumGroupsets(); gs_id++)
    {
      const auto& gs = do_problem_->GetGroupset(gs_id);
      if ((gs.first_group <= g) and (g <= gs.last_group))
      {
        groupset_ids_.push_back(gs_id);
        break;
      }
    }
  }
  else if (selected_groupset_.has_value())
  {
    const auto groupset_id = selected_groupset_.value();
    const auto& groupset = do_problem_->GetGroupsets()[groupset_id];
    for (unsigned int g = groupset.first_group; g <= groupset.last_group; ++g)
    {
      groups_.push_back(g);
      groupset_ids_.push_back(groupset_id);
    }
  }
  else
  {
    const auto& groupsets = do_problem_->GetGroupsets();
    for (unsigned int gsi = 0; gsi < groupsets.size(); ++gsi)
    {
      const auto& gs = groupsets[gsi];
      for (unsigned int g = gs.first_group; g <= gs.last_group; g++)
      {
        groups_.push_back(g);
        groupset_ids_.push_back(gsi);
      }
    }
  }
}

void
SurfacePostprocessor::Execute()
{
  for (unsigned int i = 0; i < side_faces_.size(); ++i)
  {
    switch (value_type_)
    {
      case ValueType::INTEGRAL:
        values_[i] = ComputeIntegral(side_faces_[i]);
        break;

      case ValueType::MAX:
        values_[i] = ComputeMax(side_faces_[i]);
        break;

      case ValueType::MIN:
        values_[i] = ComputeMin(side_faces_[i]);
        break;

      case ValueType::AVERAGE:
        values_[i] = ComputeAvg(side_faces_[i]);
        break;
    }
  }
}

std::vector<double>
SurfacePostprocessor::ComputeIntegral(const std::vector<SideFace>& side_faces)
{
  const auto& sdm = do_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& sweep_boundaries = do_problem_->GetSweepBoundaries();
  const auto& cell_transport_views = do_problem_->GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& groupsets = do_problem_->GetGroupsets();
  const auto& uk_man = do_problem_->GetUnknownManager();
  const auto phi = do_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  std::vector<double> local_integral(groups_.size(), 0.0);
  for (const auto [cell_id, face_idx] : side_faces)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& fe_intgrl_values = unit_cell_matrices[cell.local_id];
    const auto& IntS_shapeI = fe_intgrl_values.intS_shapeI;

    const auto& face = cell.faces[face_idx];
    const auto& bndry = sweep_boundaries.at(face.neighbor_id);
    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      const auto g = groups_[k];
      const auto& groupset = groupsets[groupset_ids_[k]];

      if (current_type_ == CurrentType::INCOMING)
        local_integral[k] +=
          ComputeInflow(groupset, cell, cell_mapping, *bndry, IntS_shapeI, face, face_idx, g);
      else if (current_type_ == CurrentType::OUTGOING)
        local_integral[k] += transport_view.GetOutflow(face_idx, g);
      else if (current_type_ == CurrentType::NET)
        local_integral[k] +=
          transport_view.GetOutflow(face_idx, g) +
          ComputeInflow(groupset, cell, cell_mapping, *bndry, IntS_shapeI, face, face_idx, g);
    }
  }

  std::vector<double> global_integral(groups_.size(), 0.0);
  for (std::size_t i = 0; i < local_integral.size(); ++i)
    mpi_comm.all_reduce(local_integral[i], global_integral[i], mpi::op::sum<double>());

  return global_integral;
}

std::vector<double>
SurfacePostprocessor::ComputeMin(const std::vector<SideFace>& side_faces)
{
  const auto& sdm = do_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& sweep_boundaries = do_problem_->GetSweepBoundaries();
  const auto& cell_transport_views = do_problem_->GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& groupsets = do_problem_->GetGroupsets();
  const auto& uk_man = do_problem_->GetUnknownManager();
  const auto phi = do_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  std::vector<double> local_values(groups_.size(), std::numeric_limits<double>::infinity());
  for (const auto [cell_id, face_idx] : side_faces)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& fe_intgrl_values = unit_cell_matrices[cell.local_id];
    const auto& IntS_shapeI = fe_intgrl_values.intS_shapeI;

    const auto& face = cell.faces[face_idx];
    const auto& bndry = sweep_boundaries.at(face.neighbor_id);
    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      const auto g = groups_[k];
      const auto& groupset = groupsets[groupset_ids_[k]];

      if (current_type_ == CurrentType::INCOMING)
        local_values[k] = std::min(
          local_values[k],
          ComputeInflow(groupset, cell, cell_mapping, *bndry, IntS_shapeI, face, face_idx, g));
      else if (current_type_ == CurrentType::OUTGOING)
        local_values[k] = std::min(local_values[k], transport_view.GetOutflow(face_idx, g));
      else if (current_type_ == CurrentType::NET)
        local_values[k] = std::min(
          local_values[k],
          transport_view.GetOutflow(face_idx, g) +
            ComputeInflow(groupset, cell, cell_mapping, *bndry, IntS_shapeI, face, face_idx, g));
    }
  }

  std::vector<double> global_values(groups_.size(), 0.0);
  for (std::size_t i = 0; i < local_values.size(); ++i)
    mpi_comm.all_reduce(local_values[i], global_values[i], mpi::op::min<double>());

  return global_values;
}

std::vector<double>
SurfacePostprocessor::ComputeMax(const std::vector<SideFace>& side_faces)
{
  const auto& sdm = do_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  const auto& sweep_boundaries = do_problem_->GetSweepBoundaries();
  const auto& cell_transport_views = do_problem_->GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& groupsets = do_problem_->GetGroupsets();
  const auto& uk_man = do_problem_->GetUnknownManager();
  const auto phi = do_problem_->GetPhiNewLocal();
  auto coord = sdm.GetSpatialWeightingFunction();

  std::vector<double> local_values(groups_.size(), -std::numeric_limits<double>::infinity());
  for (const auto [cell_id, face_idx] : side_faces)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& fe_intgrl_values = unit_cell_matrices[cell.local_id];
    const auto& IntS_shapeI = fe_intgrl_values.intS_shapeI;

    const auto& face = cell.faces[face_idx];
    const auto& bndry = sweep_boundaries.at(face.neighbor_id);
    for (std::size_t k = 0; k < groups_.size(); ++k)
    {
      const auto g = groups_[k];
      const auto& groupset = groupsets[groupset_ids_[k]];

      if (current_type_ == CurrentType::INCOMING)
        local_values[k] = std::max(
          local_values[k],
          ComputeInflow(groupset, cell, cell_mapping, *bndry, IntS_shapeI, face, face_idx, g));
      else if (current_type_ == CurrentType::OUTGOING)
        local_values[k] = std::max(local_values[k], transport_view.GetOutflow(face_idx, g));
      else if (current_type_ == CurrentType::NET)
        local_values[k] = std::max(
          local_values[k],
          transport_view.GetOutflow(face_idx, g) +
            ComputeInflow(groupset, cell, cell_mapping, *bndry, IntS_shapeI, face, face_idx, g));
    }
  }

  std::vector<double> global_values(groups_.size(), 0.0);
  for (std::size_t i = 0; i < local_values.size(); ++i)
    mpi_comm.all_reduce(local_values[i], global_values[i], mpi::op::max<double>());

  return global_values;
}

std::vector<double>
SurfacePostprocessor::ComputeAvg(const std::vector<SideFace>& side_faces)
{
  const auto& sdm = do_problem_->GetSpatialDiscretization();
  const auto& grid = sdm.GetGrid();
  auto coord = sdm.GetSpatialWeightingFunction();

  double local_weighted_area = 0.0;
  for (const auto [cell_id, face_idx] : side_faces)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto fe_face_data = cell_mapping.MakeSurfaceFiniteElementData(face_idx);

    for (const std::size_t qp : fe_face_data.GetQuadraturePointIndices())
      local_weighted_area += coord(fe_face_data.QPointXYZ(qp)) * fe_face_data.JxW(qp);
  }
  double global_weighted_area = 0.0;
  mpi_comm.all_reduce(local_weighted_area, global_weighted_area, mpi::op::sum<double>());

  auto global_weighted_integral = ComputeIntegral(side_faces);

  std::vector<double> values(groups_.size());
  for (std::size_t i = 0; i < global_weighted_integral.size(); ++i)
  {
    if (global_weighted_area > 0.0)
      values[i] = global_weighted_integral[i] / global_weighted_area;
    else
      values[i] = 0.0;
  }
  return values;
}

std::vector<std::vector<double>>
SurfacePostprocessor::GetValue() const
{
  return values_;
}

} // namespace opensn
