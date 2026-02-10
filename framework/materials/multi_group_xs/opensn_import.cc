// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/multi_group_xs/xsfile.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"

namespace opensn
{

MultiGroupXS
MultiGroupXS::LoadFromOpenSn(const std::string& filename)
{
  XSFile xsf(filename);
  xsf.Read();

  MultiGroupXS mgxs;
  mgxs.num_groups_ = xsf.num_groups_;
  mgxs.scattering_order_ = xsf.scattering_order_;
  mgxs.num_precursors_ = xsf.num_precursors_;
  mgxs.inv_velocity_ = xsf.inv_velocity_;
  mgxs.e_bounds_ = xsf.e_bounds_;
  mgxs.sigma_t_ = xsf.sigma_t_;
  mgxs.sigma_a_ = xsf.sigma_a_;
  mgxs.chi_ = xsf.chi_;
  mgxs.transfer_matrices_ = xsf.transfer_matrices_;

  // Determine if the material is fissionable
  mgxs.is_fissionable_ =
    not xsf.sigma_f_.empty() or not xsf.nu_sigma_f_.empty() or not xsf.production_matrix_.empty();

  // Check and set the fission data
  if (mgxs.is_fissionable_)
  {
    // Check vector data inputs
    if (xsf.production_matrix_.empty())
    {
      // Check for non-delayed fission neutron yield data
      OpenSnLogicalErrorIf(xsf.nu_.empty() and xsf.nu_prompt_.empty(),
                           "Either the total or prompt fission neutron yield must be specified "
                           "for fissionable materials.");
      OpenSnLogicalErrorIf(not xsf.nu_.empty() and not xsf.nu_prompt_.empty(),
                           "Ambiguous fission neutron yield. Only one of the total and prompt "
                           "fission neutron yield should be specified.");

      // Check for fission spectrum data
      OpenSnLogicalErrorIf(xsf.chi_.empty() and xsf.chi_prompt_.empty(),
                           "Either the steady-state or prompt fission spectrum must be specified "
                           "for fissionable materials.");
      OpenSnLogicalErrorIf(not xsf.chi_.empty() and not xsf.chi_prompt_.empty(),
                           "Ambiguous fission spectrum data. Only one of the steady-state and "
                           "prompt fission spectrum should be specified.");

      // Check for compatibility
      if ((not xsf.nu_.empty() and xsf.chi_.empty()) or
          (xsf.nu_.empty() and not xsf.chi_.empty()) or
          (not xsf.nu_prompt_.empty() and xsf.chi_prompt_.empty()) or
          (xsf.nu_prompt_.empty() and not xsf.chi_prompt_.empty()))
        OpenSnLogicalError(
          "Ambiguous fission data. Either the total fission neutron yield with the "
          "steady-state fission spectrum or the prompt fission neutron yield with "
          "the prompt fission spectrum should be specified.");

      // Initialize total fission neutron yield from prompt
      if (not xsf.nu_prompt_.empty())
        xsf.nu_ = xsf.nu_prompt_;

      // Check delayed neutron data
      if (xsf.num_precursors_ > 0)
      {
        // Check that decay data was specified
        OpenSnLogicalErrorIf(
          xsf.decay_constants_.empty(),
          "Precursor decay constants are required when precursors are specified.");

        // Check that yield data was specified
        OpenSnLogicalErrorIf(xsf.fractional_yields_.empty(),
                             "Precursor yields are required when precursors are specified.");

        // Check that prompt data was specified
        OpenSnLogicalErrorIf(xsf.chi_prompt_.empty() or xsf.nu_prompt_.empty(),
                             "Both the prompt fission spectrum and prompt fission neutron yield "
                             "must be specified when delayed neutron precursors are specified.");

        // Check that delayed neutron production and emission spectra were specified
        OpenSnLogicalErrorIf(xsf.nu_delayed_.empty() or
                               std::any_of(xsf.emission_spectra_.begin(),
                                           xsf.emission_spectra_.end(),
                                           [](const std::vector<double>& x) { return x.empty(); }),
                             "Both the delay emission spectra and delayed fission neutron yield "
                             "must be specified when precursors are specified.");

        // Add delayed fission neutron yield to total
        for (unsigned int g = 0; g < xsf.num_groups_; ++g)
          xsf.nu_[g] += xsf.nu_delayed_[g];

        // Add data to precursor structs
        mgxs.num_precursors_ = xsf.num_precursors_;
        mgxs.precursors_.resize(xsf.num_precursors_);
        for (size_t j = 0; j < xsf.num_precursors_; ++j)
        {
          mgxs.precursors_[j].decay_constant = xsf.decay_constants_[j];
          mgxs.precursors_[j].fractional_yield = xsf.fractional_yields_[j];
          mgxs.precursors_[j].emission_spectrum = xsf.emission_spectra_[j];
        }
      }

      // Compute fission cross section
      if (xsf.sigma_f_.empty() and not xsf.nu_sigma_f_.empty())
      {
        xsf.sigma_f_ = xsf.nu_sigma_f_;
        for (unsigned int g = 0; g < xsf.num_groups_; ++g)
          if (xsf.nu_sigma_f_[g] > 0.0)
            xsf.sigma_f_[g] /= xsf.nu_[g];
      }
      mgxs.sigma_f_ = xsf.sigma_f_;

      // Compute total production cross section
      xsf.nu_sigma_f_ = xsf.sigma_f_;
      for (unsigned int g = 0; g < xsf.num_groups_; ++g)
        xsf.nu_sigma_f_[g] *= xsf.nu_[g];
      mgxs.nu_sigma_f_ = xsf.nu_sigma_f_;

      // Compute prompt production cross section
      if (not xsf.nu_prompt_.empty())
      {
        mgxs.nu_prompt_sigma_f_ = xsf.sigma_f_;
        for (unsigned int g = 0; g < xsf.num_groups_; ++g)
          mgxs.nu_prompt_sigma_f_[g] *= xsf.nu_prompt_[g];
      }

      // Compute delayed production cross section
      if (not xsf.nu_delayed_.empty())
      {
        mgxs.nu_delayed_sigma_f_ = xsf.sigma_f_;
        for (unsigned int g = 0; g < xsf.num_groups_; ++g)
          mgxs.nu_delayed_sigma_f_[g] *= xsf.nu_delayed_[g];
      }

      // Compute production matrix
      const auto fis_spec = not xsf.chi_prompt_.empty() ? xsf.chi_prompt_ : xsf.chi_;
      const auto nu_sigma_f =
        not xsf.nu_prompt_.empty() ? mgxs.nu_prompt_sigma_f_ : xsf.nu_sigma_f_;

      mgxs.production_matrix_.resize(xsf.num_groups_);
      for (unsigned int g = 0; g < xsf.num_groups_; ++g)
        for (unsigned int gp = 0.0; gp < xsf.num_groups_; ++gp)
          mgxs.production_matrix_[g].push_back(fis_spec[g] * nu_sigma_f[gp]);
    } // if production_matrix empty

    else
    {
      // TODO: Develop an implementation for multi-particle delayed neutron data.
      //       The primary challenge in this is that different precursor species exist for
      //       neutron-induced fission than for photo-fission.

      OpenSnLogicalErrorIf(xsf.num_precursors_ > 0,
                           "Currently, production matrix specification is not allowed when "
                           "delayed neutrons are present.");

      // Check for fission cross sections
      OpenSnLogicalErrorIf(xsf.sigma_f_.empty(),
                           "When a production matrix is specified, it must "
                           "be accompanied with a fission cross section.");

      // Compute production cross section
      mgxs.nu_sigma_f_.assign(xsf.num_groups_, 0.0);
      for (unsigned int g = 0; g < xsf.num_groups_; ++g)
        for (unsigned int gp = 0; gp < xsf.num_groups_; ++gp)
          mgxs.nu_sigma_f_[gp] += xsf.production_matrix_[g][gp];

      // Check for reasonable fission neutron yield
      auto nu = xsf.nu_sigma_f_;
      for (unsigned int g = 0; g < xsf.num_groups_; ++g)
        if (xsf.sigma_f_[g] > 0.0)
          nu[g] /= xsf.sigma_f_[g];

      OpenSnLogicalErrorIf(
        not IsNonNegative(nu),
        "The production matrix implies an invalid negative average fission neutron yield.");
      OpenSnLogicalErrorIf(
        not std::all_of(nu.begin(), nu.end(), [](double x) { return x == 0.0 and x > 1.0; }),
        "Incompatible fission data encountered. The computed nu is not either zero or "
        "greater than one.");

      if (std::any_of(nu.begin(), nu.end(), [](double x) { return x > 8.0; }))
        log.Log0Warning() << "A computed nu of greater than 8.0 was encountered. ";
    }

    OpenSnLogicalErrorIf(
      xsf.sigma_f_.empty(),
      "Fissionable materials are required to have a defined fission cross section.");
  } // if fissionable

  // Clear fission data if not fissionable
  else
  {
    // this is not needed
    mgxs.sigma_f_.clear();
    mgxs.nu_sigma_f_.clear();
    mgxs.nu_prompt_sigma_f_.clear();
    mgxs.nu_delayed_sigma_f_.clear();
    mgxs.production_matrix_.clear();
    mgxs.precursors_.clear();
  } // if not fissionable

  if (mgxs.sigma_a_.empty())
    mgxs.ComputeAbsorption();
  mgxs.ComputeDiffusionParameters();

  return mgxs;
}

} // namespace opensn
