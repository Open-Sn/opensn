// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/postprocessors/cross_section_sensitivity_postprocessor.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, CrossSectionSensitivityPostprocessor);

InputParameters
CrossSectionSensitivityPostprocessor::GetInputParameters()
{
  InputParameters params;
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "A handle to an existing discrete "
                                                        "ordinates problem.");
  params.AddOptionalParameterArray(
    "block_ids",
    std::vector<int>{},
    "Block restriction for the postprocessor. Empty/unspecified means no block restriction.");
  params.AddOptionalParameterArray("logical_volumes",
                                   std::vector<std::shared_ptr<LogicalVolume>>{},
                                   "Logical volumes to restrict the computation to.");
  params.AddOptionalParameter<std::string>(
    "sensitivity_type", "sigma_t", "Sensitivity type: 'sigma_t', 'scatter', or 'production'.");
  params.AddOptionalParameter("group", 0, "Group for sigma_t or production sensitivities.");
  params.AddOptionalParameter("moment", 0, "Scattering moment/order for scatter sensitivities.");
  params.AddOptionalParameter(
    "ell", 0, "Alias for 'moment'. Specify either 'moment' or 'ell', not both.");
  params.AddOptionalParameter("from_group", 0, "Source group for scatter coefficient.");
  params.AddOptionalParameter("to_group", 0, "Destination group for scatter coefficient.");
  params.AddOptionalParameter("relative", false, "Return relative sensitivities x dR/dx.");
  params.AddOptionalParameter<std::string>(
    "forward_flux_moments", "", "Forward flux-moment file prefix. Empty uses current phi.");
  params.AddOptionalParameter<std::string>(
    "adjoint_flux_moments", "", "Adjoint flux-moment file prefix. Empty uses current phi.");
  params.AddOptionalParameter<std::string>(
    "forward_angular_fluxes", "", "Forward angular-flux file prefix. Empty uses current psi.");
  params.AddOptionalParameter<std::string>(
    "adjoint_angular_fluxes", "", "Adjoint angular-flux file prefix. Empty uses current psi.");
  params.AddOptionalParameter("flux_moments_single_file",
                              false,
                              "Read flux moments from a single file instead of one file per rank.");

  params.ConstrainParameterRange("sensitivity_type",
                                 AllowableRangeList::New({"sigma_t", "scatter", "production"}));
  return params;
}

std::shared_ptr<CrossSectionSensitivityPostprocessor>
CrossSectionSensitivityPostprocessor::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<CrossSectionSensitivityPostprocessor>(
    "lbs::CrossSectionSensitivityPostprocessor", params);
}

CrossSectionSensitivityPostprocessor::CrossSectionSensitivityPostprocessor(
  const InputParameters& params)
  : do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem")),
    block_ids_(params.GetParamVectorValue<int>("block_ids")),
    logical_volumes_(params.GetParamVectorValue<std::shared_ptr<LogicalVolume>>("logical_volumes")),
    selected_group_(params.IsParameterValid("group")
                      ? std::make_optional(params.GetParamValue<unsigned int>("group"))
                      : std::nullopt),
    ell_(params.IsParameterValid("moment")
           ? std::make_optional(params.GetParamValue<unsigned int>("moment"))
           : (params.IsParameterValid("ell")
                ? std::make_optional(params.GetParamValue<unsigned int>("ell"))
                : std::nullopt)),
    from_group_(params.IsParameterValid("from_group")
                  ? std::make_optional(params.GetParamValue<unsigned int>("from_group"))
                  : std::nullopt),
    to_group_(params.IsParameterValid("to_group")
                ? std::make_optional(params.GetParamValue<unsigned int>("to_group"))
                : std::nullopt),
    forward_flux_moments_prefix_(params.GetParamValue<std::string>("forward_flux_moments")),
    adjoint_flux_moments_prefix_(params.GetParamValue<std::string>("adjoint_flux_moments")),
    forward_angular_fluxes_prefix_(params.GetParamValue<std::string>("forward_angular_fluxes")),
    adjoint_angular_fluxes_prefix_(params.GetParamValue<std::string>("adjoint_angular_fluxes")),
    flux_moments_single_file_(params.GetParamValue<bool>("flux_moments_single_file")),
    relative_(params.GetParamValue<bool>("relative"))
{
  const auto type = params.GetParamValue<std::string>("sensitivity_type");
  if (type == "sigma_t")
    sensitivity_type_ = SensitivityType::SIGMA_T;
  else if (type == "scatter")
    sensitivity_type_ = SensitivityType::SCATTER;
  else if (type == "production")
    sensitivity_type_ = SensitivityType::PRODUCTION;
  else
    OpenSnInvalidArgument("Unsupported sensitivity_type \"" + type + "\".");

  OpenSnInvalidArgumentIf(
    selected_group_.has_value() and selected_group_.value() >= do_problem_->GetNumGroups(),
    "'group' must be less than " + std::to_string(do_problem_->GetNumGroups()) + ".");
  OpenSnInvalidArgumentIf(params.IsParameterValid("moment") and params.IsParameterValid("ell"),
                          "Specify either 'moment' or 'ell', not both.");

  CreateSpatialRestriction();
  CreateEnergyRestriction();
  CreateScatteringMomentRestriction();
  ValidateSelectedCoefficient();

  const auto num_columns =
    sensitivity_type_ == SensitivityType::SIGMA_T
      ? groups_.size()
      : (sensitivity_type_ == SensitivityType::SCATTER ? scattering_moments_.size() : 1);
  const auto n_lvs = std::max(static_cast<std::size_t>(1), logical_volumes_.size());
  values_ = NDArray<double, 2>({n_lvs, num_columns}, 0.0);
}

void
CrossSectionSensitivityPostprocessor::CreateSpatialRestriction()
{
  const auto& grid = do_problem_->GetGrid();

  if (logical_volumes_.empty())
  {
    cell_local_ids_.resize(1);
    for (const auto& cell : grid->local_cells)
      if (block_ids_.empty() or
          std::find(block_ids_.begin(), block_ids_.end(), cell.block_id) != block_ids_.end())
        cell_local_ids_[0].push_back(cell.local_id);
  }
  else
  {
    cell_local_ids_.resize(logical_volumes_.size());
    for (unsigned int i = 0; i < logical_volumes_.size(); ++i)
      cell_local_ids_[i] = GetLogicalVolumeCellIDs(logical_volumes_[i]);
  }

  OpenSnLogicalErrorIf(cell_local_ids_.empty(),
                       "No mesh cells were selected by a CrossSectionSensitivityPostprocessor.");
}

std::vector<std::uint32_t>
CrossSectionSensitivityPostprocessor::GetLogicalVolumeCellIDs(
  std::shared_ptr<LogicalVolume> log_vol) const
{
  const auto& grid = do_problem_->GetGrid();
  std::vector<std::uint32_t> cell_ids;
  for (const auto& cell : grid->local_cells)
  {
    if (not log_vol->Inside(cell.centroid))
      continue;
    if (block_ids_.empty() or
        std::find(block_ids_.begin(), block_ids_.end(), cell.block_id) != block_ids_.end())
      cell_ids.push_back(cell.local_id);
  }
  return cell_ids;
}

void
CrossSectionSensitivityPostprocessor::CreateEnergyRestriction()
{
  if (selected_group_.has_value())
  {
    groups_.push_back(selected_group_.value());
  }
  else
  {
    for (unsigned int g = 0; g < do_problem_->GetNumGroups(); ++g)
      groups_.push_back(g);
  }
}

void
CrossSectionSensitivityPostprocessor::CreateScatteringMomentRestriction()
{
  if (sensitivity_type_ != SensitivityType::SCATTER)
    return;

  if (ell_.has_value())
  {
    scattering_moments_.push_back(ell_.value());
    return;
  }

  for (unsigned int ell = 0; ell <= do_problem_->GetScatteringOrder(); ++ell)
    scattering_moments_.push_back(ell);
}

void
CrossSectionSensitivityPostprocessor::ValidateSelectedCoefficient() const
{
  if (sensitivity_type_ == SensitivityType::SIGMA_T)
    return;

  if (sensitivity_type_ == SensitivityType::PRODUCTION)
  {
    OpenSnInvalidArgumentIf(not selected_group_.has_value(),
                            "'group' is required for production sensitivities.");
    OpenSnInvalidArgumentIf(selected_group_.value() >= do_problem_->GetNumGroups(),
                            "'group' must be a valid group index.");
    return;
  }

  OpenSnInvalidArgumentIf(not from_group_.has_value() or not to_group_.has_value(),
                          "'from_group' and 'to_group' are required for scatter sensitivities.");
  OpenSnInvalidArgumentIf(from_group_.value() >= do_problem_->GetNumGroups() or
                            to_group_.value() >= do_problem_->GetNumGroups(),
                          "'from_group' and 'to_group' must be valid group indices.");

  OpenSnInvalidArgumentIf(scattering_moments_.empty(), "No scattering moments were selected.");
  for (const auto ell : scattering_moments_)
    OpenSnInvalidArgumentIf(ell >= do_problem_->GetNumMoments(),
                            "Selected scattering moment must be less than the number of flux "
                            "moments.");
}

void
CrossSectionSensitivityPostprocessor::LoadAngularFluxes(
  std::vector<std::vector<double>>& forward_psi, std::vector<std::vector<double>>& adjoint_psi)
{
  if (forward_angular_fluxes_prefix_.empty())
    forward_psi = do_problem_->GetPsiNewLocal();
  else
    DiscreteOrdinatesProblemIO::ReadAngularFluxes(
      *do_problem_, forward_angular_fluxes_prefix_, forward_psi);

  if (adjoint_angular_fluxes_prefix_.empty())
    adjoint_psi = do_problem_->GetPsiNewLocal();
  else
    DiscreteOrdinatesProblemIO::ReadAngularFluxes(
      *do_problem_, adjoint_angular_fluxes_prefix_, adjoint_psi);

  OpenSnInvalidArgumentIf(forward_psi.empty() or adjoint_psi.empty(),
                          "sigma_t sensitivities require forward and adjoint angular fluxes.");
  OpenSnInvalidArgumentIf(forward_psi.size() != adjoint_psi.size(),
                          "Forward and adjoint angular fluxes have different groupset counts.");
  for (size_t gs = 0; gs < forward_psi.size(); ++gs)
    OpenSnInvalidArgumentIf(forward_psi[gs].size() != adjoint_psi[gs].size(),
                            "Forward and adjoint angular flux groupset " + std::to_string(gs) +
                              " have different sizes.");
}

void
CrossSectionSensitivityPostprocessor::LoadFluxMoments(std::vector<double>& forward_phi,
                                                      std::vector<double>& adjoint_phi)
{
  if (forward_flux_moments_prefix_.empty())
    forward_phi = do_problem_->GetPhiNewLocal();
  else
    LBSSolverIO::ReadFluxMoments(
      *do_problem_, forward_flux_moments_prefix_, flux_moments_single_file_, forward_phi);

  if (adjoint_flux_moments_prefix_.empty())
    adjoint_phi = do_problem_->GetPhiNewLocal();
  else
    LBSSolverIO::ReadFluxMoments(
      *do_problem_, adjoint_flux_moments_prefix_, flux_moments_single_file_, adjoint_phi);

  OpenSnInvalidArgumentIf(forward_phi.empty() or adjoint_phi.empty(),
                          "scatter and production sensitivities require forward and adjoint flux "
                          "moments.");
  OpenSnInvalidArgumentIf(forward_phi.size() != adjoint_phi.size(),
                          "Forward and adjoint flux moments have different sizes.");
}

void
CrossSectionSensitivityPostprocessor::Execute()
{
  if (sensitivity_type_ == SensitivityType::SIGMA_T)
  {
    std::vector<std::vector<double>> forward_psi;
    std::vector<std::vector<double>> adjoint_psi;
    LoadAngularFluxes(forward_psi, adjoint_psi);

    for (size_t i = 0; i < cell_local_ids_.size(); ++i)
    {
      const auto sensitivities =
        ComputeTotalSensitivity(cell_local_ids_[i], forward_psi, adjoint_psi);
      for (size_t j = 0; j < sensitivities.size(); ++j)
        values_(i, j) = sensitivities[j];
    }
  }
  else
  {
    std::vector<double> forward_phi;
    std::vector<double> adjoint_phi;
    LoadFluxMoments(forward_phi, adjoint_phi);

    for (size_t i = 0; i < cell_local_ids_.size(); ++i)
    {
      const auto sensitivities =
        sensitivity_type_ == SensitivityType::SCATTER
          ? ComputeScatterSensitivity(cell_local_ids_[i], forward_phi, adjoint_phi)
          : ComputeProductionSensitivity(cell_local_ids_[i], forward_phi, adjoint_phi);
      for (size_t j = 0; j < sensitivities.size(); ++j)
        values_(i, j) = sensitivities[j];
    }
  }
}

void
CrossSectionSensitivityPostprocessor::ApplyKEigenvalueScaling(const double k_eff)
{
  OpenSnInvalidArgumentIf(k_eff <= 0.0, "The forward k-eigenvalue must be positive.");

  std::vector<double> forward_phi;
  std::vector<double> adjoint_phi;
  LoadFluxMoments(forward_phi, adjoint_phi);

  ScaleForKEigenvalueSensitivity(k_eff, ComputeFissionDenominator(forward_phi, adjoint_phi));
}

std::vector<double>
CrossSectionSensitivityPostprocessor::ComputeTotalSensitivity(
  const std::vector<std::uint32_t>& cell_local_ids,
  const std::vector<std::vector<double>>& forward_psi,
  const std::vector<std::vector<double>>& adjoint_psi) const
{
  const auto& grid = do_problem_->GetGrid();
  const auto& discretization = do_problem_->GetSpatialDiscretization();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& groupsets = do_problem_->GetGroupsets();

  std::vector<int> group_to_column(do_problem_->GetNumGroups(), -1);
  for (size_t k = 0; k < groups_.size(); ++k)
    group_to_column[groups_[k]] = static_cast<int>(k);

  std::vector<double> local(groups_.size(), 0.0);
  for (const auto cell_id : cell_local_ids)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);
    const auto& sigma_t = xs->GetSigmaTotal();

    for (const auto& groupset : groupsets)
    {
      const auto gs = groupset.id;
      const auto& uk_man = groupset.psi_uk_man_;
      const auto& quadrature = groupset.quadrature;
      const auto num_gs_groups = groupset.GetNumGroups();
      const auto num_gs_angles = quadrature->GetOmegas().size();

      for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
      {
        const auto V_i = fe_values.intV_shapeI(i);
        for (size_t n = 0; n < num_gs_angles; ++n)
        {
          const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, n, 0);
          const auto weight = quadrature->GetWeight(n) * V_i;
          for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
          {
            const auto g = groupset.first_group + gsg;
            const int column = group_to_column[g];
            if (column < 0)
              continue;
            double contribution =
              -weight * adjoint_psi[gs][dof_map + gsg] * forward_psi[gs][dof_map + gsg];
            if (relative_)
              contribution *= sigma_t[g];
            local[column] += contribution;
          }
        }
      }
    }
  }

  std::vector<double> global(groups_.size(), 0.0);
  for (size_t k = 0; k < local.size(); ++k)
    mpi_comm.all_reduce(local[k], global[k], mpi::op::sum<double>());
  return global;
}

std::vector<double>
CrossSectionSensitivityPostprocessor::ComputeScatterSensitivity(
  const std::vector<std::uint32_t>& cell_local_ids,
  const std::vector<double>& forward_phi,
  const std::vector<double>& adjoint_phi) const
{
  const auto& grid = do_problem_->GetGrid();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& transport_views = do_problem_->GetCellTransportViews();
  const auto from_group = from_group_.value_or(0);
  const auto to_group = to_group_.value_or(0);

  const LBSGroupset* coefficient_groupset = nullptr;
  for (const auto& groupset : do_problem_->GetGroupsets())
    if (to_group >= groupset.first_group and to_group <= groupset.last_group)
    {
      coefficient_groupset = &groupset;
      break;
    }
  OpenSnLogicalErrorIf(coefficient_groupset == nullptr,
                       "Unable to identify groupset for scattering destination group.");
  const auto& moment_map = coefficient_groupset->quadrature->GetMomentToHarmonicsIndexMap();

  std::vector<int> ell_to_column(do_problem_->GetNumMoments(), -1);
  for (size_t k = 0; k < scattering_moments_.size(); ++k)
    ell_to_column[scattering_moments_[k]] = static_cast<int>(k);

  std::vector<double> local(scattering_moments_.size(), 0.0);
  for (const auto cell_id : cell_local_ids)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = transport_views[cell.local_id];
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);

    for (unsigned int m = 0; m < moment_map.size(); ++m)
    {
      const auto ell = moment_map[m].ell;
      const int column = ell < ell_to_column.size() ? ell_to_column[ell] : -1;
      if (column < 0)
        continue;

      double relative_coefficient = 1.0;
      if (relative_)
      {
        if (ell < xs->GetTransferMatrices().size())
          relative_coefficient = xs->GetTransferMatrix(ell).GetValueIJ(to_group, from_group);
        else
          relative_coefficient = 0.0;
      }

      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const auto dof_map = transport_view.MapDOF(i, m, 0);
        local[column] += relative_coefficient * adjoint_phi[dof_map + to_group] *
                         forward_phi[dof_map + from_group] * fe_values.intV_shapeI(i);
      }
    }
  }

  std::vector<double> global(scattering_moments_.size(), 0.0);
  for (size_t k = 0; k < local.size(); ++k)
    mpi_comm.all_reduce(local[k], global[k], mpi::op::sum<double>());
  return global;
}

std::vector<double>
CrossSectionSensitivityPostprocessor::ComputeProductionSensitivity(
  const std::vector<std::uint32_t>& cell_local_ids,
  const std::vector<double>& forward_phi,
  const std::vector<double>& adjoint_phi) const
{
  const auto& grid = do_problem_->GetGrid();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& transport_views = do_problem_->GetCellTransportViews();
  const auto group = selected_group_.value_or(0);

  double local = 0.0;
  for (const auto cell_id : cell_local_ids)
  {
    const auto& cell = grid->local_cells[cell_id];
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = transport_views[cell.local_id];
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);

    if (not xs->IsFissionable())
      continue;

    const auto& chi = xs->GetChi();
    const auto& nu_sigma_f = xs->GetNuSigmaF();

    OpenSnLogicalErrorIf(group >= nu_sigma_f.size(),
                         "Production sensitivity group is outside the fission XS data.");
    OpenSnLogicalErrorIf(chi.size() < do_problem_->GetNumGroups(),
                         "The fission spectrum does not contain all energy groups.");

    // Sensitivity is with respect to nu_sigma_f (see post_processors.rst: absolute =
    // dR/d(nu_sigma_f), relative = nu_sigma_f * dR/d(nu_sigma_f)), not sigma_f - the
    // loop below already accumulates dR/d(nu_sigma_f) directly, so the absolute case
    // needs no extra factor. Scaling by nu_sigma_f/sigma_f (i.e. nu, holding nu fixed)
    // would instead give dR/d(sigma_f), a different coefficient than the one this
    // function is documented to report.
    const double coefficient = relative_ ? nu_sigma_f[group] : 1.0;

    if (coefficient == 0.0)
      continue;

    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const auto dof_map = transport_view.MapDOF(i, 0, 0);

      double adjoint_fission_importance = 0.0;
      for (unsigned int g = 0; g < do_problem_->GetNumGroups(); ++g)
        adjoint_fission_importance += chi[g] * adjoint_phi[dof_map + g];

      local += coefficient * forward_phi[dof_map + group] * adjoint_fission_importance *
               fe_values.intV_shapeI(i);
    }
  }

  double global = 0.0;
  mpi_comm.all_reduce(local, global, mpi::op::sum<double>());
  return {global};
}

double
CrossSectionSensitivityPostprocessor::ComputeFissionDenominator(
  const std::vector<double>& forward_phi, const std::vector<double>& adjoint_phi) const
{
  const auto& grid = do_problem_->GetGrid();
  const auto& unit_cell_matrices = do_problem_->GetUnitCellMatrices();
  const auto& transport_views = do_problem_->GetCellTransportViews();

  double local = 0.0;
  for (const auto& cell : grid->local_cells)
  {
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);
    if (not xs->IsFissionable())
      continue;

    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = transport_views[cell.local_id];
    const auto& chi = xs->GetChi();
    const auto& nu_sigma_f = xs->GetNuSigmaF();

    OpenSnLogicalErrorIf(chi.size() < do_problem_->GetNumGroups() or
                           nu_sigma_f.size() < do_problem_->GetNumGroups(),
                         "Fission XS data does not contain all energy groups.");

    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const auto dof_map = transport_view.MapDOF(i, 0, 0);

      double adjoint_fission_importance = 0.0;
      double forward_fission_production = 0.0;
      for (unsigned int g = 0; g < do_problem_->GetNumGroups(); ++g)
      {
        adjoint_fission_importance += chi[g] * adjoint_phi[dof_map + g];
        forward_fission_production += nu_sigma_f[g] * forward_phi[dof_map + g];
      }

      local += adjoint_fission_importance * forward_fission_production * fe_values.intV_shapeI(i);
    }
  }

  double global = 0.0;
  mpi_comm.all_reduce(local, global, mpi::op::sum<double>());
  return global;
}

void
CrossSectionSensitivityPostprocessor::ScaleForKEigenvalueSensitivity(const double k_eff,
                                                                     const double denominator)
{
  OpenSnInvalidArgumentIf(std::abs(denominator) == 0.0,
                          "The fission normalization denominator <psi_adj,F psi> is zero.");

  const double scale = sensitivity_type_ == SensitivityType::PRODUCTION
                         ? k_eff / denominator
                         : k_eff * k_eff / denominator;
  const auto num_columns =
    sensitivity_type_ == SensitivityType::SIGMA_T
      ? groups_.size()
      : (sensitivity_type_ == SensitivityType::SCATTER ? scattering_moments_.size()
                                                       : static_cast<std::size_t>(1));

  for (size_t i = 0; i < cell_local_ids_.size(); ++i)
    for (size_t j = 0; j < num_columns; ++j)
      values_(i, j) *= scale;
}

const NDArray<double, 2>&
CrossSectionSensitivityPostprocessor::GetValue() const
{
  return values_;
}

} // namespace opensn
