// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/hdf_utils.h"
#include "framework/utils/utils.h"
#include <numeric>
#include <algorithm>
#include <cmath>

namespace opensn
{

namespace
{

bool
IsFinite(const std::vector<double>& values)
{
  return std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); });
}

/// Reads a one-dimensional dataset of `size` finite, non-negative values.
std::vector<double>
ReadOpenMCVector(hid_t file_id,
                 const std::string& path,
                 const std::string& name,
                 const size_t size,
                 const std::string& file_name)
{
  std::vector<double> values;
  OpenSnLogicalErrorIf(not H5ReadDataset1D<double>(file_id, path + name, values),
                       "Failed to read the one-dimensional dataset \"" + path + name + "\" from " +
                         file_name + ".");
  OpenSnLogicalErrorIf(values.size() != size,
                       "Dataset \"" + path + name + "\" in " + file_name + " has " +
                         std::to_string(values.size()) + " values; expected " +
                         std::to_string(size) + ".");
  OpenSnLogicalErrorIf(not IsFinite(values) or not IsNonNegative(values),
                       "Dataset \"" + path + name + "\" in " + file_name +
                         " must be finite and non-negative.");
  return values;
}

/**
 * Reads a [delayed group][group] dataset of finite, non-negative values. With one delayed group,
 * a one-dimensional [group] dataset is accepted too.
 */
std::vector<std::vector<double>>
ReadOpenMCFamilyData(hid_t file_id,
                     const std::string& path,
                     const std::string& name,
                     const unsigned int num_families,
                     const unsigned int num_groups,
                     const std::string& file_name)
{
  std::vector<std::vector<double>> values;
  if (not H5ReadDataset2D<double>(file_id, path + name, values))
  {
    std::vector<double> values_1d;
    OpenSnLogicalErrorIf(not H5ReadDataset1D<double>(file_id, path + name, values_1d),
                         "Failed to read dataset \"" + path + name + "\" from " + file_name +
                           ". Delayed-neutron data must be a [delayed group][group] array "
                           "(matrix and angular representations are not supported).");
    OpenSnLogicalErrorIf(num_families != 1,
                         "Dataset \"" + path + name + "\" in " + file_name +
                           " is one-dimensional, but the file has " + std::to_string(num_families) +
                           " delayed groups.");
    values.push_back(std::move(values_1d));
  }
  OpenSnLogicalErrorIf(values.size() != num_families or values.front().size() != num_groups,
                       "Dataset \"" + path + name + "\" in " + file_name + " must have shape [" +
                         std::to_string(num_families) + "][" + std::to_string(num_groups) + "].");
  const bool valid = std::all_of(values.begin(),
                                 values.end(),
                                 [](const std::vector<double>& family)
                                 { return IsFinite(family) and IsNonNegative(family); });
  OpenSnLogicalErrorIf(not valid,
                       "Dataset \"" + path + name + "\" in " + file_name +
                         " must be finite and non-negative.");
  return values;
}

/// Normalizes a non-negative spectrum to unit sum.
void
NormalizeOpenMCSpectrum(std::vector<double>& spectrum, const std::string& description)
{
  const double sum = std::accumulate(spectrum.begin(), spectrum.end(), 0.0);
  OpenSnLogicalErrorIf(not(sum > 0.0), description + " must have at least one non-zero value.");
  for (auto& value : spectrum)
    value /= sum;
}

} // namespace

MultiGroupXS
MultiGroupXS::LoadFromOpenMC(const std::string& file_name,
                             const std::string& dataset_name,
                             double temperature,
                             const std::vector<std::string>& extra_xs_names)
{
  MultiGroupXS mgxs;

  // Open file
  const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (file.Id() < 0)
  {
    std::string err_msg = "Unable to open " + file_name + " or it is not a valid HDF5 file.\n";
    throw std::logic_error(err_msg);
  }

  // Check file type
  std::string filetype;
  H5ReadAttribute<std::string>(file.Id(), "filetype", filetype);
  if (filetype != "mgxs")
    throw std::runtime_error(file_name + " is not a valid OpenMC library file");

  log.Log() << "Reading OpenMC cross-section file \"" << file_name << "\"\n";

  // Number of groups
  if (not H5ReadAttribute<unsigned int>(file.Id(), "energy_groups", mgxs.num_groups_))
    throw std::runtime_error("Failure reading \"energy_groups\" from " + file_name);

  // Group structure
  H5ReadDataset1D<double>(file.Id(), "/group structure", mgxs.e_bounds_);
  std::reverse(mgxs.e_bounds_.begin(), mgxs.e_bounds_.end());

  // Temperature
  log.Log0() << file_name + " cross-section data evaluated at " << temperature << "K\n";
  mgxs.temperature_ = temperature;

  // Base path
  std::stringstream ss;
  ss << std::fixed << std::setprecision(0) << mgxs.temperature_ << "K";
  std::string path = "/" + dataset_name + "/" + ss.str() + "/";
  if (!H5Has(file.Id(), path))
  {
    throw std::runtime_error("Could not find dataset " + dataset_name + "/" + ss.str() + " in " +
                             file_name);
  }
  mgxs.temperature_ = temperature;

  // Scattering order
  if (not H5ReadGroupAttribute<unsigned int>(
        file.Id(), dataset_name, "order", mgxs.scattering_order_))
    throw std::runtime_error("Failure reading \"order\" from " + file_name);

  // Scattering shape
  std::string scatter_shape;
  H5ReadGroupAttribute<std::string>(file.Id(), dataset_name, "scatter_shape", scatter_shape);
  if (scatter_shape != "[G][G'][Order]")
    throw std::runtime_error(file_name + " has an unsupported scatter shape");

  // Number of precursors
  int delayed_groups = 0;
  OpenSnLogicalErrorIf(
    not H5ReadOptionalAttribute<int>(file.Id(), "delayed_groups", delayed_groups),
    "Failure reading \"delayed_groups\" from " + file_name + ".");
  OpenSnLogicalErrorIf(delayed_groups < 0,
                       "\"delayed_groups\" in " + file_name + " must be non-negative.");
  mgxs.num_precursors_ = static_cast<unsigned int>(delayed_groups);

  // Inverse velocity
  H5ReadDataset1D<double>(file.Id(), path + "inverse-velocity", mgxs.inv_velocity_);
  OpenSnLogicalErrorIf(not IsNonNegative(mgxs.inv_velocity_),
                       "Only positive inverse velocity values are permitted.");

  // Total
  H5ReadDataset1D<double>(file.Id(), path + "total", mgxs.sigma_t_);
  OpenSnLogicalErrorIf(mgxs.sigma_t_.empty(),
                       "\"total\" data block not found in " + file_name + ".");
  OpenSnLogicalErrorIf(not IsNonNegative(mgxs.sigma_t_),
                       "Only non-negative total cross-section values are permitted.");

  // Absorption
  // We do not read OpenMC absorption values (even when available) but prefer to re-compute them
  // using the total XS and the transfer matrix XS that were just read in order to avoid issues due
  // to small statistical noise

  // Transfer
  // Note that the scatter matrices in OpenMC are stored as the transposed version of what OpenSn
  // uses.
  if (H5Has(file.Id(), path + "scatter_data/scatter_matrix"))
  {
    mgxs.transfer_matrices_.assign(mgxs.scattering_order_ + 1,
                                   SparseMatrix(mgxs.num_groups_, mgxs.num_groups_));
    std::vector<double> flat_scatter_matrix;
    H5ReadDataset1D<double>(file.Id(), path + "scatter_data/scatter_matrix", flat_scatter_matrix);
    std::vector<int> g_min, g_max;
    H5ReadDataset1D<int>(file.Id(), path + "scatter_data/g_min", g_min);
    H5ReadDataset1D<int>(file.Id(), path + "scatter_data/g_max", g_max);
    int fidx = 0;
    for (unsigned int gp = 0; gp < mgxs.num_groups_; ++gp)
      for (int g = g_min[gp]; g <= g_max[gp]; ++g)
        for (unsigned int n = 0; n < mgxs.scattering_order_ + 1; ++n, ++fidx)
          mgxs.transfer_matrices_.at(n).Insert(g - 1, gp, flat_scatter_matrix[fidx]);
  }

  if (mgxs.sigma_a_.empty())
    mgxs.ComputeAbsorption();
  mgxs.ComputeDiffusionParameters();

  // Is fissionable?
  H5ReadGroupAttribute<bool>(file.Id(), dataset_name, "fissionable", mgxs.is_fissionable_);

  // OpenMC gives the delayed-neutron data of a library with delayed groups in one of two vector
  // forms. As in OpenMC, the beta form is used when "beta" exists.
  //  - beta: "nu-fission" and "chi" with "beta", the delayed fraction of each delayed group, either
  //    [delayed group] or [delayed group][group]. Prompt and delayed neutrons share chi.
  //  - prompt/delayed: "prompt-nu-fission", "delayed-nu-fission" [delayed group][group],
  //    "chi-prompt", and "chi-delayed" [delayed group][group]. "nu-fission" is their sum.
  const unsigned int num_families = mgxs.num_precursors_;
  const unsigned int num_groups = mgxs.num_groups_;
  const bool has_delayed_data = mgxs.is_fissionable_ and num_families > 0;
  const bool beta_form = has_delayed_data and H5Has(file.Id(), path + "beta");
  std::vector<std::vector<double>> delayed_nu_sigma_f_by_family;
  std::vector<double> nu_prompt_sigma_f;

  if (mgxs.is_fissionable_)
  {
    // Fission
    mgxs.sigma_f_ = ReadOpenMCVector(file.Id(), path, "fission", num_groups, file_name);
    if (not HasNonZero(mgxs.sigma_f_))
    {
      log.Log0Warning() << "The fission cross section specified in "
                        << "\"" << file_name << "\" is uniformly zero... Clearing it.";
      mgxs.sigma_f_.clear();
    }

    // Nu-Fission
    if (has_delayed_data and not beta_form)
    {
      delayed_nu_sigma_f_by_family = ReadOpenMCFamilyData(
        file.Id(), path, "delayed-nu-fission", num_families, num_groups, file_name);
      std::vector<double> nu_delayed_sigma_f(num_groups, 0.0);
      for (const auto& family : delayed_nu_sigma_f_by_family)
        for (unsigned int g = 0; g < num_groups; ++g)
          nu_delayed_sigma_f[g] += family[g];

      if (H5Has(file.Id(), path + "prompt-nu-fission"))
      {
        nu_prompt_sigma_f =
          ReadOpenMCVector(file.Id(), path, "prompt-nu-fission", num_groups, file_name);
        mgxs.nu_sigma_f_.assign(num_groups, 0.0);
        for (unsigned int g = 0; g < num_groups; ++g)
          mgxs.nu_sigma_f_[g] = nu_prompt_sigma_f[g] + nu_delayed_sigma_f[g];
      }
      else
      {
        // Without prompt-nu-fission, the prompt part is the rest of the total.
        mgxs.nu_sigma_f_ = ReadOpenMCVector(file.Id(), path, "nu-fission", num_groups, file_name);
        nu_prompt_sigma_f = mgxs.nu_sigma_f_;
        for (unsigned int g = 0; g < num_groups; ++g)
          nu_prompt_sigma_f[g] -= nu_delayed_sigma_f[g];
        OpenSnLogicalErrorIf(not IsNonNegative(nu_prompt_sigma_f),
                             "\"delayed-nu-fission\" in " + file_name +
                               " exceeds \"nu-fission\" in some group.");
      }
    }
    else
    {
      mgxs.nu_sigma_f_ = ReadOpenMCVector(file.Id(), path, "nu-fission", num_groups, file_name);
    }
    OpenSnLogicalErrorIf(
      not IsNonNegative(mgxs.nu_sigma_f_),
      "Only non-negative total fission multiplication cross-section values are permitted.");
    if (not HasNonZero(mgxs.nu_sigma_f_))
    {
      log.Log0Warning() << "The production cross section specified in "
                        << "\"" << file_name << "\" is uniformly zero... Clearing it.";
      mgxs.nu_sigma_f_.clear();
    }

    // A material flagged as fissionable but with uniformly-zero fission and production
    // cross sections has no actual fission capability.
    if (mgxs.sigma_f_.empty() or mgxs.nu_sigma_f_.empty())
      mgxs.is_fissionable_ = false;
  }

  if (mgxs.is_fissionable_)
  {
    // Chi, the fission spectrum. The prompt/delayed form does not need it.
    if (has_delayed_data and not beta_form and not H5Has(file.Id(), path + "chi"))
      mgxs.chi_.clear();
    else
    {
      H5ReadDataset1D<double>(file.Id(), path + "chi", mgxs.chi_);
      OpenSnLogicalErrorIf(mgxs.chi_.empty(), "\"chi\" data block not found in " + file_name + ".");
      OpenSnLogicalErrorIf(mgxs.chi_.size() != num_groups,
                           "\"chi\" in " + file_name + " must have one value per group.");
      OpenSnLogicalErrorIf(not IsNonNegative(mgxs.chi_),
                           "Steady-state fission spectrum must be non-negative.");
      NormalizeOpenMCSpectrum(mgxs.chi_, "Steady-state fission spectrum");
    }

    // Nu (computed)
    auto nu = mgxs.nu_sigma_f_;
    for (unsigned int g = 0; g < num_groups; ++g)
      if (mgxs.sigma_f_[g] > 0.0)
        nu[g] /= mgxs.sigma_f_[g];
    OpenSnLogicalErrorIf(not IsNonNegative(nu),
                         "The fission data imply an invalid negative average fission neutron "
                         "yield.");
    OpenSnLogicalErrorIf(
      not std::all_of(nu.begin(), nu.end(), [](double x) { return x == 0.0 or x > 1.0; }),
      "Incompatible fission data encountered. The computed nu is not either zero or "
      "greater than one.");
    if (std::any_of(nu.begin(), nu.end(), [](double x) { return x > 8.0; }))
      log.Log0Warning() << "A computed nu of greater than 8.0 was encountered. ";

    if (has_delayed_data)
    {
      std::vector<double> chi_prompt;
      std::vector<std::vector<double>> chi_delayed_by_family;
      if (beta_form)
      {
        // Delayed fraction of each delayed group: [delayed group][group], or [delayed group]
        // when it does not depend on the fissioning group
        std::vector<std::vector<double>> beta;
        if (H5ReadDataset2D<double>(file.Id(), path + "beta", beta))
          beta = ReadOpenMCFamilyData(file.Id(), path, "beta", num_families, num_groups, file_name);
        else
        {
          const auto beta_by_family =
            ReadOpenMCVector(file.Id(), path, "beta", num_families, file_name);
          for (const double beta_j : beta_by_family)
            beta.emplace_back(num_groups, beta_j);
        }

        delayed_nu_sigma_f_by_family.assign(num_families, std::vector<double>(num_groups, 0.0));
        nu_prompt_sigma_f = mgxs.nu_sigma_f_;
        for (unsigned int j = 0; j < num_families; ++j)
          for (unsigned int g = 0; g < num_groups; ++g)
          {
            delayed_nu_sigma_f_by_family[j][g] = beta[j][g] * mgxs.nu_sigma_f_[g];
            nu_prompt_sigma_f[g] -= delayed_nu_sigma_f_by_family[j][g];
          }
        OpenSnLogicalErrorIf(not IsNonNegative(nu_prompt_sigma_f),
                             "\"beta\" in " + file_name +
                               " sums to more than one over the delayed groups.");
        chi_prompt = mgxs.chi_;
        chi_delayed_by_family.assign(num_families, mgxs.chi_);
      }
      else
      {
        chi_delayed_by_family =
          ReadOpenMCFamilyData(file.Id(), path, "chi-delayed", num_families, num_groups, file_name);
        for (auto& family : chi_delayed_by_family)
          NormalizeOpenMCSpectrum(family, "Each delayed group's \"chi-delayed\"");

        if (H5Has(file.Id(), path + "chi-prompt"))
          chi_prompt = ReadOpenMCVector(file.Id(), path, "chi-prompt", num_groups, file_name);
        else
        {
          OpenSnLogicalErrorIf(mgxs.chi_.empty(),
                               "Neither \"chi-prompt\" nor \"chi\" found in " + file_name + ".");
          log.Log0Warning() << "OpenMC delayed-neutron data found in \"" << file_name
                            << "\" without chi-prompt. Using chi as the prompt spectrum.";
          chi_prompt = mgxs.chi_;
        }
        NormalizeOpenMCSpectrum(chi_prompt, "\"chi-prompt\"");
      }

      auto decay_constants =
        ReadOpenMCVector(file.Id(), path, "decay-rate", num_families, file_name);
      OpenSnLogicalErrorIf(not IsPositive(decay_constants),
                           "\"decay-rate\" in " + file_name + " must be positive.");

      // OpenSn gives each precursor family one fractional yield, its share of the delayed
      // nu-fission summed over the groups. This reproduces the data only when the share does not
      // depend on the fissioning group.
      std::vector<double> nu_delayed_sigma_f(num_groups, 0.0);
      std::vector<double> family_totals(num_families, 0.0);
      for (unsigned int j = 0; j < num_families; ++j)
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          nu_delayed_sigma_f[g] += delayed_nu_sigma_f_by_family[j][g];
          family_totals[j] += delayed_nu_sigma_f_by_family[j][g];
        }
      const double delayed_total = std::accumulate(family_totals.begin(), family_totals.end(), 0.0);

      double max_share_deviation = 0.0;
      mgxs.precursors_.resize(num_families);
      for (unsigned int j = 0; j < num_families; ++j)
      {
        auto& precursor = mgxs.precursors_[j];
        precursor.decay_constant = decay_constants[j];
        precursor.emission_spectrum = chi_delayed_by_family[j];
        precursor.fractional_yield = delayed_total > 0.0 ? family_totals[j] / delayed_total : 0.0;
        for (unsigned int g = 0; g < num_groups; ++g)
          if (nu_delayed_sigma_f[g] > 0.0)
            max_share_deviation =
              std::max(max_share_deviation,
                       std::fabs(delayed_nu_sigma_f_by_family[j][g] / nu_delayed_sigma_f[g] -
                                 precursor.fractional_yield));
      }
      if (max_share_deviation > 1.0e-6)
        log.Log0Warning() << "The delayed groups' shares of the delayed nu-fission in \""
                          << file_name << "\" (" << dataset_name
                          << ") depend on the fissioning group (by up to " << max_share_deviation
                          << "). OpenSn uses one fractional yield per precursor family, its "
                             "share of the delayed nu-fission summed over the groups, so the "
                             "delayed neutron production and precursor concentrations are "
                             "approximate.";

      mgxs.nu_prompt_sigma_f_ = std::move(nu_prompt_sigma_f);
      mgxs.nu_delayed_sigma_f_ = std::move(nu_delayed_sigma_f);

      // The production matrix is prompt only; the precursors give the delayed production.
      mgxs.production_matrix_.assign(num_groups, std::vector<double>(num_groups, 0.0));
      for (unsigned int gp = 0; gp < num_groups; ++gp)
        for (unsigned int g = 0; g < num_groups; ++g)
          mgxs.production_matrix_[g][gp] = chi_prompt[g] * mgxs.nu_prompt_sigma_f_[gp];
    }
    else
    {
      // Production matrix (computed), total production
      mgxs.production_matrix_.resize(num_groups, std::vector<double>(num_groups));
      for (unsigned int gp = 0; gp < num_groups; ++gp)
        for (unsigned int g = 0; g < num_groups; ++g)
          mgxs.production_matrix_[g][gp] = mgxs.chi_[g] * mgxs.nu_sigma_f_[gp];
    }
  } // if fissionable
  else
  {
    // Clear fission data if not fissionable
    mgxs.num_precursors_ = 0;
    mgxs.sigma_f_.clear();
    mgxs.nu_sigma_f_.clear();
    mgxs.nu_prompt_sigma_f_.clear();
    mgxs.nu_delayed_sigma_f_.clear();
    mgxs.production_matrix_.clear();
    mgxs.precursors_.clear();
  }

  for (const auto& xs_name : extra_xs_names)
  {
    if (!H5Has(file.Id(), path + xs_name))
    {
      std::string msg = "Requested XS \"";
      msg += xs_name;
      msg += "\" not found in ";
      msg += file_name;
      throw std::runtime_error(msg);
    }
    std::vector<double> xs_vals;
    if (!H5ReadDataset1D<double>(file.Id(), path + xs_name, xs_vals))
    {
      std::string msg = "Failed to read XS \"";
      msg += xs_name;
      msg += "\" from ";
      msg += file_name;
      throw std::runtime_error(msg);
    }
    if (xs_vals.size() != mgxs.num_groups_)
      throw std::runtime_error("Requested XS \"" + xs_name + "\" does not have num_groups entries");
    mgxs.custom_xs_[xs_name] = std::move(xs_vals);
  }

  return mgxs;
}

} // namespace opensn
