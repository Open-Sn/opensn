// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/physics/physics_material/mgxs/mgxs.h"
#include "framework/logging/log.h"

namespace opensn
{

void
MGXS::Reset()
{
  num_groups_ = 0;
  scattering_order_ = 0;
  num_precursors_ = 0;
  is_fissionable_ = false;

  sigma_t_.clear();
  sigma_a_.clear();
  transfer_matrices_.clear();

  sigma_f_.clear();
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

  // Monte-Carlo quantities
  cdf_gprime_g_.clear();
  scat_angles_gprime_g_.clear();
}

void
MGXS::ComputeAbsorption()
{
  sigma_a_.assign(num_groups_, 0.0);

  // Compute for a pure absorber
  if (transfer_matrices_.empty())
  {
    for (size_t g = 0; g < num_groups_; ++g)
      sigma_a_[g] = sigma_t_[g];
  }

  // Estimate from a transfer matrix
  else
  {
    log.Log0Warning() << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices_.front();
    for (size_t g = 0; g < num_groups_; ++g)
    {
      // Estimate the scattering cross section
      double sigma_s = 0.0;
      for (size_t row = 0; row < S0.NumRows(); ++row)
      {
        const auto& cols = S0.rowI_indices_[row];
        const auto& vals = S0.rowI_values_[row];
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
  }   // if scattering present
}

void
MGXS::ComputeDiffusionParameters()
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
  for (size_t g = 0; g < num_groups_; ++g)
  {
    // Determine transport correction
    double sigma_1 = 0.0;
    if (S.size() > 1)
    {
      for (size_t gp = 0; gp < num_groups_; ++gp)
      {
        const auto& cols = S[1].rowI_indices_[gp];
        const auto& vals = S[1].rowI_values_[gp];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sigma_1 += vals[t];
            break;
          }
      } // for gp
    }   // if moment 1 available

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
    diffusion_coeff_[g] = std::fmin(1.0e12, 1.0 / 3.0 / (sigma_t_[g] - sigma_1));

    // Determine within group scattering
    if (not S.empty())
    {
      const auto& cols = S[0].rowI_indices_[g];
      const auto& vals = S[0].rowI_values_[g];
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
MGXS::SetScalingFactor(const double factor)
{
  const double m = factor / scaling_factor_;
  scaling_factor_ = factor;

  // Apply to STL vector-based data
  for (size_t g = 0; g < num_groups_; ++g)
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
    for (size_t g = 0; g < num_groups_; ++g)
      for (const auto& [_, gp, sig_ell] : S_ell.Row(g))
        sig_ell *= m;

  // Reinitialize diffusion
  diffusion_initialized_ = false;
  ComputeDiffusionParameters();
}

void
MGXS::TransposeTransferAndProduction()
{
  // Transpose transfer matrices
  for (unsigned int ell = 0; ell <= scattering_order_; ++ell)
  {
    const auto& S_ell = transfer_matrices_[ell];
    SparseMatrix S_ell_transpose(num_groups_, num_groups_);
    for (size_t g = 0; g < num_groups_; ++g)
    {
      const size_t row_len = S_ell.rowI_indices_[g].size();
      const size_t* col_ptr = S_ell.rowI_indices_[g].data();
      const double* val_ptr = S_ell.rowI_values_[g].data();

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
    for (size_t g = 0; g < num_groups_; ++g)
      for (size_t gp = 0; gp < num_groups_; ++gp)
        transposed_production_matrix_[g].push_back(F[gp][g]);
  }
}

} // namespace opensn
