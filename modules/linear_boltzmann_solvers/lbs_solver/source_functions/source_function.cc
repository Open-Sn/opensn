// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/source_functions/source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "caliper/cali.h"

namespace opensn
{

SourceFunction::SourceFunction(const LBSSolver& lbs_solver) : lbs_solver_(lbs_solver)
{
}

void
SourceFunction::operator()(const LBSGroupset& groupset,
                           std::vector<double>& q,
                           const std::vector<double>& phi,
                           const SourceFlags source_flags)
{
  CALI_CXX_MARK_SCOPE("SourceFunction::operator");

  if (source_flags.Empty())
    return;

  apply_fixed_src_ = (source_flags & APPLY_FIXED_SOURCES);
  apply_wgs_scatter_src_ = (source_flags & APPLY_WGS_SCATTER_SOURCES);
  apply_ags_scatter_src_ = (source_flags & APPLY_AGS_SCATTER_SOURCES);
  apply_wgs_fission_src_ = (source_flags & APPLY_WGS_FISSION_SOURCES);
  apply_ags_fission_src_ = (source_flags & APPLY_AGS_FISSION_SOURCES);
  suppress_wg_scatter_src_ = (source_flags & SUPPRESS_WG_SCATTER);

  const auto& densities = lbs_solver_.GetDensitiesLocal();

  // Get group setup
  gs_i_ = static_cast<size_t>(groupset.groups.front().id);
  gs_f_ = static_cast<size_t>(groupset.groups.back().id);

  first_grp_ = static_cast<size_t>(lbs_solver_.GetGroups().front().id);
  last_grp_ = static_cast<size_t>(lbs_solver_.GetGroups().back().id);

  default_zero_src_.assign(lbs_solver_.GetGroups().size(), 0.0);

  const auto& cell_transport_views = lbs_solver_.GetCellTransportViews();

  const auto num_moments = lbs_solver_.GetNumMoments();
  const auto& ext_src_moments_local = lbs_solver_.GetExtSrcMomentsLocal();

  const auto& m_to_ell_em_map = groupset.quadrature->GetMomentToHarmonicsIndexMap();

  // Apply all nodal sources
  const auto& grid = lbs_solver_.GetGrid();
  for (const auto& cell : grid.local_cells)
  {
    const auto& rho = densities[cell.local_id];
    const auto& transport_view = cell_transport_views[cell.local_id];
    cell_volume_ = transport_view.GetVolume();

    // Obtain xs
    const auto& xs = transport_view.GetXS();

    const auto& S = xs.GetTransferMatrices();
    const auto& F = xs.GetProductionMatrix();
    const auto& precursors = xs.GetPrecursors();
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    // Loop over nodes
    const auto num_nodes = transport_view.GetNumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      // Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments); ++m)
      {
        const auto ell = m_to_ell_em_map[m].ell;
        const auto uk_map = transport_view.MapDOF(i, m, 0);
        const double* phi_im = &phi[uk_map];

        // Declare moment src
        fixed_src_moments_ = default_zero_src_.data();

        if (lbs_solver_.GetOptions().use_src_moments)
          fixed_src_moments_ = &ext_src_moments_local[uk_map];

        // Loop over groupset groups
        for (size_t g = gs_i_; g <= gs_f_; ++g)
        {
          g_ = g;

          double rhs = 0.0;

          // Apply fixed sources
          if (apply_fixed_src_)
            rhs += this->AddSourceMoments();

          // Apply scattering sources
          if (ell < S.size())
          {
            const auto& S_ell = S[ell];
            // Add Across GroupSet Scattering (AGS)
            if (apply_ags_scatter_src_)
              for (const auto& [_, gp, sigma_sm] : S_ell.Row(g))
                if (gp < gs_i_ or gp > gs_f_)
                  rhs += rho * sigma_sm * phi_im[gp];

            // Add Within GroupSet Scattering (WGS)
            if (apply_wgs_scatter_src_)
              for (const auto& [_, gp, sigma_sm] : S_ell.Row(g))
                if (gp >= gs_i_ and gp <= gs_f_)
                {
                  if (suppress_wg_scatter_src_ and g_ == gp)
                    continue;
                  rhs += rho * sigma_sm * phi_im[gp];
                }
          }

          // Apply fission sources
          if (xs.IsFissionable() and ell == 0)
          {
            const auto& F_g = F[g];
            if (apply_ags_fission_src_)
              for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
                if (gp < gs_i_ or gp > gs_f_)
                  rhs += rho * F_g[gp] * phi_im[gp];

            if (apply_wgs_fission_src_)
              for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
                rhs += rho * F_g[gp] * phi_im[gp];

            if (lbs_solver_.GetOptions().use_precursors)
              rhs += this->AddDelayedFission(precursors, rho, nu_delayed_sigma_f, &phi[uk_map]);
          }

          // Add to destination vector
          q[uk_map + g] += rhs;

        } // for g
      }   // for m
    }     // for dof i
  }       // for cell

  AddAdditionalSources(groupset, q, phi, source_flags);
}

double
SourceFunction::AddSourceMoments() const
{
  return fixed_src_moments_[g_];
}

double
SourceFunction::AddDelayedFission(const PrecursorList& precursors,
                                  const double& rho,
                                  const std::vector<double>& nu_delayed_sigma_f,
                                  const double* phi) const
{
  double value = 0.0;
  if (apply_ags_fission_src_)
    for (size_t gp = first_grp_; gp <= last_grp_; ++gp)
      if (gp < gs_i_ or gp > gs_f_)
        for (const auto& precursor : precursors)
          value += precursor.emission_spectrum[g_] * precursor.fractional_yield * rho *
                   nu_delayed_sigma_f[gp] * phi[gp];

  if (apply_wgs_fission_src_)
    for (size_t gp = gs_i_; gp <= gs_f_; ++gp)
      for (const auto& precursor : precursors)
        value += precursor.emission_spectrum[g_] * precursor.fractional_yield * rho *
                 nu_delayed_sigma_f[gp] * phi[gp];

  return value;
}

void
SourceFunction::AddPointSources(const LBSGroupset& groupset,
                                std::vector<double>& q,
                                const std::vector<double>&,
                                const SourceFlags source_flags)
{
  const bool apply_fixed_src = (source_flags & APPLY_FIXED_SOURCES);

  const auto& transport_views = lbs_solver_.GetCellTransportViews();

  const auto gs_i = groupset.groups.front().id;
  const auto gs_f = groupset.groups.back().id;

  // Apply point sources
  if (not lbs_solver_.GetOptions().use_src_moments and apply_fixed_src)
  {
    for (const auto& point_source : lbs_solver_.GetPointSources())
    {
      for (const auto& subscriber : point_source->GetSubscribers())
      {
        auto& transport_view = transport_views[subscriber.cell_local_id];

        const auto& strength = point_source->GetStrength();
        const auto& node_weights = subscriber.node_weights;
        const auto volume_weight = subscriber.volume_weight;

        for (size_t i = 0; i < transport_view.GetNumNodes(); ++i)
        {
          const auto uk_map = transport_view.MapDOF(i, 0, 0);
          for (size_t g = gs_i; g <= gs_f; ++g)
            q[uk_map + g] += strength[g] * node_weights(i) * volume_weight;
        } // for node i
      }   // for subscriber
    }     // for point source
  }
}

void
SourceFunction::AddVolumetricSources(const LBSGroupset& groupset,
                                     std::vector<double>& q,
                                     const std::vector<double>& phi,
                                     const SourceFlags source_flags)
{
  const bool apply_fixed_src = source_flags & APPLY_FIXED_SOURCES;

  const auto& grid = lbs_solver_.GetGrid();
  const auto& discretization = lbs_solver_.GetSpatialDiscretization();
  const auto& cell_transport_views = lbs_solver_.GetCellTransportViews();
  const auto num_groups = lbs_solver_.GetNumGroups();

  const auto gs_i = groupset.groups.front().id;
  const auto gs_f = groupset.groups.back().id;

  // Go through each volumetric source, and its subscribing cells
  if (not lbs_solver_.GetOptions().use_src_moments and apply_fixed_src)
  {
    for (const auto& volumetric_source : lbs_solver_.GetVolumetricSources())
    {
      for (const auto local_id : volumetric_source->GetSubscribers())
      {
        const auto& cell = grid.local_cells[local_id];
        const auto& transport_view = cell_transport_views[local_id];
        const auto nodes = discretization.GetCellNodeLocations(cell);
        const auto num_cell_nodes = discretization.GetCellNumNodes(cell);

        // Go through each of the cell nodes
        for (size_t i = 0; i < num_cell_nodes; ++i)
        {
          // Compute group-wise values for this node
          const auto src = (*volumetric_source)(cell, nodes[i], num_groups);

          // Contribute to the source moments
          const auto dof_map = transport_view.MapDOF(i, 0, 0);
          for (size_t g = gs_i; g <= gs_f; ++g)
            q[dof_map + g] += src[g];
        } // for node i
      }   // for subscriber
    }     // for volumetric source
  }
}

} // namespace opensn
