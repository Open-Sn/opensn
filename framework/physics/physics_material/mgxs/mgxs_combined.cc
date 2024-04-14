// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/physics/physics_material/mgxs/mgxs.h"
#include "framework/logging/log.h"

namespace opensn
{

void
MGXS::Initialize(std::vector<std::pair<int, double>>& combinations)
{
  Reset();

  // Pickup all xs and make sure they are valid
  std::vector<std::shared_ptr<MGXS>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;
  double Nf_total = 0.0; // Total density of fissile materials

  // Loop over cross sections
  for (auto combo : combinations)
  {
    // Get the cross section from the stack
    std::shared_ptr<MGXS> xs;
    xs = GetStackItemPtr(multigroup_xs_stack, combo.first, std::string(__FUNCTION__));
    xsecs.push_back(xs);

    // Set the scaling factor
    SetScalingFactor(combo.second);

    // Increment densities
    if (xs->IsFissionable())
    {
      is_fissionable_ = true;
      Nf_total += xs->ScalingFactor();
    }

    // Define and check number of groups
    if (xsecs.size() == 1)
      n_grps = xs->NumGroups();
    OpenSnLogicalErrorIf(xs->NumGroups() != n_grps,
                         "All cross sections being combined must have the same group structure.");

    // Increment number of precursors
    n_precs += num_precursors_;
  } // for cross section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.
  if (Nf_total < 1.0e-12)
    is_fissionable_ = false;

  // Check to ensure that all fissionable cross sections contain either
  // prompt/delayed fission data or total fission data
  if (n_precs > 0)
    for (const auto& xs : xsecs)
      OpenSnLogicalErrorIf(xs->IsFissionable() and xs->NumPrecursors() == 0,
                           "If precursors are specified, all fissionable cross sections must "
                           "specify precursors.");

  // Initialize the data
  num_groups_ = n_grps;
  scattering_order_ = 0;
  num_precursors_ = n_precs;
  for (const auto& xs : xsecs)
    scattering_order_ = std::max(scattering_order_, xs->ScatteringOrder());

  // Init mandatory cross sections
  sigma_t_.assign(n_grps, 0.0);
  sigma_a_.assign(n_grps, 0.0);

  // Init transfer matrices only if at least one exists
  if (std::any_of(xsecs.begin(),
                  xsecs.end(),
                  [](const std::shared_ptr<MGXS>& x) { return not x->TransferMatrices().empty(); }))
    transfer_matrices_.assign(scattering_order_ + 1, SparseMatrix(num_groups_, num_groups_));

  // Init fission data
  if (is_fissionable_)
  {
    sigma_f_.assign(n_grps, 0.0);
    nu_sigma_f_.assign(n_grps, 0.0);
    production_matrix_.assign(num_groups_, std::vector<double>(num_groups_, 0.0));

    // Init prompt/delayed fission data
    if (n_precs > 0)
    {
      nu_prompt_sigma_f_.assign(n_grps, 0.0);
      nu_delayed_sigma_f_.assign(n_grps, 0.0);
      precursors_.resize(n_precs);
    }
  }

  // Combine the data
  size_t precursor_count = 0;
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    // Fraction of fissile density
    const auto N_i = xsecs[x]->ScalingFactor();
    const auto ff_i = xsecs[x]->IsFissionable() ? N_i / Nf_total : 0.0;

    // Combine cross sections
    const auto& sig_t = xsecs[x]->SigmaTotal();
    const auto& sig_a = xsecs[x]->SigmaAbsorption();
    const auto& sig_f = xsecs[x]->SigmaFission();
    const auto& nu_p_sig_f = xsecs[x]->NuPromptSigmaF();
    const auto& nu_d_sig_f = xsecs[x]->NuDelayedSigmaF();
    const auto& F = xsecs[x]->ProductionMatrix();

    // Here, raw cross sections are scaled by densities and spectra by
    // fractional densities. The latter is done to preserve a unit spectra.
    for (size_t g = 0; g < n_grps; ++g)
    {
      sigma_t_[g] += sig_t[g];
      sigma_a_[g] += sig_a[g];

      if (xsecs[x]->IsFissionable())
      {
        sigma_f_[g] += sig_f[g];
        nu_sigma_f_[g] += sig_f[g];
        for (size_t gp = 0; gp < num_groups_; ++gp)
          production_matrix_[g][gp] += F[g][gp];

        if (n_precs > 0)
        {
          nu_prompt_sigma_f_[g] += nu_p_sig_f[g];
          nu_delayed_sigma_f_[g] += nu_d_sig_f[g];
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

    if (xsecs[x]->NumPrecursors() > 0)
    {
      const auto& precursors = xsecs[x]->Precursors();
      for (size_t j = 0; j < xsecs[x]->NumPrecursors(); ++j)
      {
        size_t count = precursor_count + j;
        const auto& precursor = precursors[j];
        precursors_[count].decay_constant = precursor.decay_constant;
        precursors_[count].fractional_yield = precursor.fractional_yield * ff_i;
        precursors_[count].emission_spectrum = precursor.emission_spectrum;
      } // for j

      precursor_count += xsecs[x]->NumPrecursors();
    }

    // Set inverse velocity data
    if (x == 0 and xsecs[x]->InverseVelocity().empty())
      inv_velocity_ = xsecs[x]->InverseVelocity();
    OpenSnLogicalErrorIf(
      xsecs[x]->InverseVelocity() != inv_velocity_,
      "All cross sections being combined must have the same group-wise velocities.");

    // Combine transfer matrices

    // This step is somewhat tricky. The cross sections aren't guaranteed
    // to have the same sparsity patterns and therefore simply adding them
    // together has to take the sparse matrix's protection mechanisms into
    // account.

    if (not xsecs[x]->TransferMatrices().empty())
    {
      for (size_t m = 0; m < xsecs[x]->ScatteringOrder() + 1; ++m)
      {
        auto& Sm = transfer_matrices_[m];
        const auto& Sm_other = xsecs[x]->TransferMatrix(m);
        for (size_t g = 0; g < num_groups_; ++g)
        {
          const auto& cols = Sm_other.rowI_indices_[g];
          const auto& vals = Sm_other.rowI_values_[g];
          for (size_t t = 0; t < cols.size(); ++t)
            Sm.InsertAdd(g, t, vals[t]);
        }
      }
    }
  } // for cross sections

  ComputeDiffusionParameters();
}

} // namespace opensn
