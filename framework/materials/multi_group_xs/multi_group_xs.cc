// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

MultiGroupXS
MultiGroupXS::CreateSimpleOneGroup(double sigma_t, double c, double velocity)
{
  MultiGroupXS mgxs;

  mgxs.num_groups_ = 1;
  mgxs.sigma_t_.resize(mgxs.num_groups_, sigma_t);
  mgxs.sigma_a_.resize(mgxs.num_groups_, sigma_t * (1.0 - c));
  mgxs.transfer_matrices_.emplace_back(mgxs.num_groups_, mgxs.num_groups_);
  auto& S = mgxs.transfer_matrices_.back();
  S.SetDiagonal(std::vector<double>(mgxs.num_groups_, sigma_t * c));

  if (velocity > 0.0)
    mgxs.inv_velocity_.assign(mgxs.num_groups_, 1.0 / velocity);

  mgxs.ComputeDiffusionParameters();

  return mgxs;
}

MultiGroupXS
MultiGroupXS::Combine(
  const std::vector<std::pair<std::shared_ptr<MultiGroupXS>, double>>& combinations)
{
  MultiGroupXS mgxs;

  // Pickup all xs and make sure they are valid
  std::vector<std::shared_ptr<MultiGroupXS>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;
  double Nf_total = 0.0; // Total density of fissile materials

  // Loop over cross sections
  for (const auto& combo : combinations)
  {
    // Get the cross section
    auto xs = combo.first;
    xsecs.push_back(xs);

    // Set the scaling factor
    xs->SetScalingFactor(combo.second);

    // Increment densities
    if (xs->IsFissionable())
    {
      mgxs.is_fissionable_ = true;
      Nf_total += xs->GetScalingFactor();
    }

    // Define and check number of groups
    if (xsecs.size() == 1)
      n_grps = xs->GetNumGroups();
    OpenSnLogicalErrorIf(xs->GetNumGroups() != n_grps,
                         "All cross sections being combined must have the same group structure.");

    // Increment number of precursors
    n_precs += xs->GetNumPrecursors();
  } // for cross section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.
  if (Nf_total < 1.0e-12)
    mgxs.is_fissionable_ = false;

  // Check to ensure that all fissionable cross sections contain either
  // prompt/delayed fission data or total fission data
  if (n_precs > 0)
    for (const auto& xs : xsecs)
      OpenSnLogicalErrorIf(xs->IsFissionable() and xs->GetNumPrecursors() == 0,
                           "If precursors are specified, all fissionable cross sections must "
                           "specify precursors.");

  // Initialize the data
  mgxs.num_groups_ = n_grps;
  mgxs.scattering_order_ = 0;
  mgxs.num_precursors_ = n_precs;
  for (const auto& xs : xsecs)
    mgxs.scattering_order_ = std::max(mgxs.scattering_order_, xs->GetScatteringOrder());

  // Init mandatory cross sections
  mgxs.sigma_t_.assign(n_grps, 0.0);
  mgxs.sigma_a_.assign(n_grps, 0.0);

  // Init transfer matrices only if at least one exists
  if (std::any_of(xsecs.begin(),
                  xsecs.end(),
                  [](const std::shared_ptr<MultiGroupXS>& x)
                  { return not x->GetTransferMatrices().empty(); }))
    mgxs.transfer_matrices_.assign(mgxs.scattering_order_ + 1,
                                   SparseMatrix(mgxs.num_groups_, mgxs.num_groups_));

  // Init fission data
  if (mgxs.is_fissionable_)
  {
    mgxs.sigma_f_.assign(n_grps, 0.0);
    mgxs.chi_.assign(n_grps, 0.0);
    mgxs.nu_sigma_f_.assign(n_grps, 0.0);
    mgxs.production_matrix_.assign(mgxs.num_groups_, std::vector<double>(mgxs.num_groups_, 0.0));

    // Init prompt/delayed fission data
    if (n_precs > 0)
    {
      mgxs.nu_prompt_sigma_f_.assign(n_grps, 0.0);
      mgxs.nu_delayed_sigma_f_.assign(n_grps, 0.0);
      mgxs.precursors_.resize(n_precs);
    }
  }

  // Combine the data
  size_t precursor_count = 0;
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    // Fraction of fissile density
    const auto N_i = xsecs[x]->GetScalingFactor();
    const auto ff_i = xsecs[x]->IsFissionable() ? N_i / Nf_total : 0.0;

    // Combine cross sections
    const auto& sig_t = xsecs[x]->GetSigmaTotal();
    const auto& sig_a = xsecs[x]->GetSigmaAbsorption();
    const auto& chi = xsecs[x]->GetChi();
    const auto& sig_f = xsecs[x]->GetSigmaFission();
    const auto& nu_p_sig_f = xsecs[x]->GetNuPromptSigmaF();
    const auto& nu_d_sig_f = xsecs[x]->GetNuDelayedSigmaF();
    const auto& F = xsecs[x]->GetProductionMatrix();

    // Here, raw cross sections are scaled by densities and spectra by
    // fractional densities. The latter is done to preserve a unit spectra.
    for (unsigned int g = 0; g < n_grps; ++g)
    {
      mgxs.sigma_t_[g] += sig_t[g];
      mgxs.sigma_a_[g] += sig_a[g];

      if (xsecs[x]->IsFissionable())
      {
        mgxs.sigma_f_[g] += sig_f[g];
        if (not chi.empty())
          mgxs.chi_[g] += ff_i * chi[g];
        mgxs.nu_sigma_f_[g] += sig_f[g];
        for (unsigned int gp = 0; gp < mgxs.num_groups_; ++gp)
          mgxs.production_matrix_[g][gp] += F[g][gp];

        if (n_precs > 0)
        {
          mgxs.nu_prompt_sigma_f_[g] += nu_p_sig_f[g];
          mgxs.nu_delayed_sigma_f_[g] += nu_d_sig_f[g];
        }
      }
    } // for g

    // Combine precursor data
    // Here, all precursors across all materials are stored. The decay constants and delayed
    // spectrum are what they are, however, some special treatment must be given to the yields.
    // Because the yield dictates what fraction of delayed neutrons are produced from a
    // given family, the sum over all families must yield unity. To achieve this end, all
    // precursor yields must be scaled based on the fraction of the total density of materials
    // with precursors they make up.

    if (xsecs[x]->GetNumPrecursors() > 0)
    {
      const auto& precursors = xsecs[x]->GetPrecursors();
      for (size_t j = 0; j < xsecs[x]->GetNumPrecursors(); ++j)
      {
        size_t count = precursor_count + j;
        const auto& precursor = precursors[j];
        mgxs.precursors_[count].decay_constant = precursor.decay_constant;
        mgxs.precursors_[count].fractional_yield = precursor.fractional_yield * ff_i;
        mgxs.precursors_[count].emission_spectrum = precursor.emission_spectrum;
      } // for j

      precursor_count += xsecs[x]->GetNumPrecursors();
    }

    // Set inverse velocity data
    if (x == 0 and not xsecs[x]->GetInverseVelocity().empty())
      mgxs.inv_velocity_ = xsecs[x]->GetInverseVelocity();
    if (not mgxs.inv_velocity_.empty())
      OpenSnLogicalErrorIf(
        xsecs[x]->GetInverseVelocity() != mgxs.inv_velocity_,
        "All cross sections being combined must have the same group-wise velocities.");

    // Combine transfer matrices

    // This step is somewhat tricky. The cross sections aren't guaranteed
    // to have the same sparsity patterns and therefore simply adding them
    // together has to take the sparse matrix's protection mechanisms into
    // account.

    if (not xsecs[x]->GetTransferMatrices().empty())
    {
      for (unsigned int m = 0; m < xsecs[x]->GetScatteringOrder() + 1; ++m)
      {
        auto& Sm = mgxs.transfer_matrices_[m];
        const auto& Sm_other = xsecs[x]->GetTransferMatrix(m);
        for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
        {
          const auto& cols = Sm_other.rowI_indices[g];
          const auto& vals = Sm_other.rowI_values[g];
          for (size_t t = 0; t < cols.size(); ++t)
            Sm.InsertAdd(g, t, vals[t]);
        }
      }
    }
  } // for cross sections

  mgxs.ComputeDiffusionParameters();

  return mgxs;
}

void
MultiGroupXS::Reset()
{
  num_groups_ = 0;
  scattering_order_ = 0;
  num_precursors_ = 0;
  is_fissionable_ = false;

  sigma_t_.clear();
  sigma_a_.clear();
  transfer_matrices_.clear();

  sigma_f_.clear();
  chi_.clear();
  nu_sigma_f_.clear();
  nu_prompt_sigma_f_.clear();
  nu_delayed_sigma_f_.clear();
  production_matrix_.clear();
  precursors_.clear();

  inv_velocity_.clear();

  // Diffusion quantities
  diffusion_initialized_ = false;
  sigma_tr_.clear();
  diffusion_coeff_.clear();
  sigma_r_.clear();
  sigma_s_gtog_.clear();
  custom_xs_.clear();
}

bool
MultiGroupXS::HasCustomXS(const std::string& name) const
{
  return custom_xs_.find(name) != custom_xs_.end();
}

const std::vector<double>&
MultiGroupXS::GetCustomXS(const std::string& name) const
{
  const auto it = custom_xs_.find(name);
  if (it == custom_xs_.end())
    throw std::runtime_error("MultiGroupXS: Custom XS not found: " + name);
  return it->second;
}

std::vector<std::string>
MultiGroupXS::GetCustomXSNames() const
{
  std::vector<std::string> names;
  names.reserve(custom_xs_.size());
  for (const auto& entry : custom_xs_)
    names.push_back(entry.first);
  return names;
}

void
MultiGroupXS::ComputeAbsorption()
{
  sigma_a_.assign(num_groups_, 0.0);

  // Compute for a pure absorber
  if (transfer_matrices_.empty())
  {
    for (unsigned int g = 0; g < num_groups_; ++g)
      sigma_a_[g] = sigma_t_[g];
  }

  // Estimate from a transfer matrix
  else
  {
    log.Log0Warning() << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices_.front();
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      // Estimate the scattering cross section
      double sigma_s = 0.0;
      for (size_t row = 0; row < S0.GetNumRows(); ++row)
      {
        const auto& cols = S0.rowI_indices[row];
        const auto& vals = S0.rowI_values[row];
        for (size_t t = 0; t < cols.size(); ++t)
        {
          if (cols[t] == g)
          {
            sigma_s += vals[t];
            break;
          }
        }
      }
      sigma_a_[g] = sigma_t_[g] - sigma_s;

      // TODO: Should negative absorption be allowed?
      if (sigma_a_[g] < 0.0)
        log.Log0Warning() << "Negative absorption cross section encountered in group " << g
                          << " when estimating from the transfer matrices";
    } // for g
  } // if scattering present
}

void
MultiGroupXS::ComputeDiffusionParameters()
{
  if (diffusion_initialized_)
    return;

  // Initialize diffusion data
  sigma_tr_.resize(num_groups_, 0.0);
  diffusion_coeff_.resize(num_groups_, 1.0);
  sigma_s_gtog_.resize(num_groups_, 0.0);
  sigma_r_.resize(num_groups_, 0.0);

  // Perform computations group-wise
  const auto& S = transfer_matrices_;
  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    // Determine transport correction
    double sigma_1 = 0.0;
    if (S.size() > 1)
    {
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
      {
        const auto& cols = S[1].rowI_indices[gp];
        const auto& vals = S[1].rowI_values[gp];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sigma_1 += vals[t];
            break;
          }
      } // for gp
    } // if moment 1 available

    // Compute transport cross section
    if (sigma_1 >= sigma_t_[g])
    {
      log.Log0Warning() << "Negative transport cross section found for "
                        << "group " << g << " in call to " << __FUNCTION__ << ". "
                        << "sigma_t=" << sigma_t_[g] << " sigma_1=" << sigma_1 << ". "
                        << "Setting sigma_1=0, sigma_tr=sigma_t for this group.";
      sigma_1 = 0.0;
    }
    sigma_tr_[g] = sigma_t_[g] - sigma_1;

    // Compute the diffusion coefficient
    // Cap the value for when sigma_t - sigma_1 is near zero
    if ((sigma_t_[g] - sigma_1) < 1.0e-12)
      diffusion_coeff_[g] = 1.0e12;
    else
      diffusion_coeff_[g] = 1.0 / 3.0 / (sigma_t_[g] - sigma_1);

    // Determine within group scattering
    if (not S.empty())
    {
      const auto& cols = S[0].rowI_indices[g];
      const auto& vals = S[0].rowI_values[g];
      for (size_t t = 0; t < cols.size(); ++t)
      {
        if (cols[t] == g)
        {
          sigma_s_gtog_[g] = vals[t];
          break;
        }
      }
    }

    // Compute removal cross section
    sigma_r_[g] = std::max(0.0, sigma_t_[g] - sigma_s_gtog_[g]);
  } // for g

  diffusion_initialized_ = true;
}

void
MultiGroupXS::SetScalingFactor(const double factor)
{
  const double m = factor / scaling_factor_;
  scaling_factor_ = factor;

  // Apply to STL vector-based data
  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    sigma_t_[g] *= m;
    sigma_a_[g] *= m;

    if (is_fissionable_)
    {
      sigma_f_[g] *= m;
      nu_sigma_f_[g] *= m;
      if (num_precursors_ > 0)
      {
        nu_prompt_sigma_f_[g] *= m;
        nu_delayed_sigma_f_[g] *= m;
      }

      for (auto& x : production_matrix_[g])
        x *= m;
    }
  }

  // Apply to transfer matrices
  for (auto& S_ell : transfer_matrices_)
    for (unsigned int g = 0; g < num_groups_; ++g)
      for (const auto& [_, gp, sig_ell] : S_ell.Row(g))
        sig_ell *= m;

  // Reinitialize diffusion
  diffusion_initialized_ = false;
  ComputeDiffusionParameters();
}

void
MultiGroupXS::TransposeTransferAndProduction()
{
  // Transpose transfer matrices
  for (unsigned int ell = 0; ell <= scattering_order_; ++ell)
  {
    const auto& S_ell = transfer_matrices_[ell];
    SparseMatrix S_ell_transpose(num_groups_, num_groups_);
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      const size_t row_len = S_ell.rowI_indices[g].size();
      const size_t* col_ptr = S_ell.rowI_indices[g].data();
      const double* val_ptr = S_ell.rowI_values[g].data();

      for (size_t j = 0; j < row_len; ++j)
        S_ell_transpose.Insert(*col_ptr++, g, *val_ptr++);
    }
    transposed_transfer_matrices_.push_back(S_ell_transpose);
  } // for ell

  // Transpose production matrices
  if (is_fissionable_)
  {
    transposed_production_matrix_.clear();
    transposed_production_matrix_.resize(num_groups_);
    const auto& F = production_matrix_;
    for (unsigned int g = 0; g < num_groups_; ++g)
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
        transposed_production_matrix_[g].push_back(F[gp][g]);
  }
}

} // namespace opensn
