// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/hdf_utils.h"
#include "framework/utils/utils.h"
#include <numeric>
#include <algorithm>
#include <cstdint>

namespace opensn
{

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
  unsigned int delayed_groups = 0;
  if (not H5ReadAttribute<unsigned int>(file.Id(), "delayed_groups", delayed_groups))
    delayed_groups = 0;
  OpenSnLogicalErrorIf(delayed_groups < 0, "OpenMC number of delayed_groups must be non-negative.");
  mgxs.num_precursors_ = delayed_groups;

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
  if (mgxs.is_fissionable_)
  {
    // Fission
    H5ReadDataset1D<double>(file.Id(), path + "fission", mgxs.sigma_f_);
    OpenSnLogicalErrorIf(mgxs.sigma_f_.empty(),
                         "\"fission\" data block not found in " + file_name + ".");
    OpenSnLogicalErrorIf(not IsNonNegative(mgxs.sigma_f_),
                         "Only non-negative fission cross-section values are permitted.");
    if (not HasNonZero(mgxs.sigma_f_))
    {
      log.Log0Warning() << "The fission cross section specified in "
                        << "\"" << file_name << "\" is uniformly zero... Clearing it.";
      mgxs.sigma_f_.clear();
    }

    // Nu-Fission
    H5ReadDataset1D<double>(file.Id(), path + "nu-fission", mgxs.nu_sigma_f_);
    OpenSnLogicalErrorIf(mgxs.nu_sigma_f_.empty(),
                         "\"nu-fission\" data block not found in " + file_name + ".");
    OpenSnLogicalErrorIf(
      not IsNonNegative(mgxs.nu_sigma_f_),
      "Only non-negative total fission multiplication cross-section values are permitted.");
    if (not HasNonZero(mgxs.nu_sigma_f_))
    {
      log.Log0Warning() << "The production cross section specified in "
                        << "\"" << file_name << "\" is uniformly zero... Clearing it.";
      mgxs.nu_sigma_f_.clear();
    }

    // Chi
    H5ReadDataset1D<double>(file.Id(), path + "chi", mgxs.chi_);
    OpenSnLogicalErrorIf(mgxs.chi_.empty(), "\"chi\" data block not found in " + file_name + ".");
    OpenSnLogicalErrorIf(not HasNonZero(mgxs.chi_),
                         "Steady-state fission spectrum must have at least one non-zero value.");
    OpenSnLogicalErrorIf(not IsNonNegative(mgxs.chi_),
                         "Steady-state fission spectrum must be non-negative.");
    // Normalizing
    const auto sum = std::accumulate(mgxs.chi_.begin(), mgxs.chi_.end(), 0.0);
    std::transform(
      mgxs.chi_.begin(), mgxs.chi_.end(), mgxs.chi_.begin(), [sum](double& x) { return x / sum; });

    // Nu (computed)
    auto nu = mgxs.nu_sigma_f_;
    for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
      if (mgxs.sigma_f_[g] > 0.0)
        nu[g] /= mgxs.sigma_f_[g];

    if (mgxs.num_precursors_ > 0)
    {
      std::vector<double> nu_prompt_sigma_f;
      if (H5Has(file.Id(), path + "prompt-nu-fission"))
        H5ReadDataset1D<double>(file.Id(), path + "prompt-nu-fission", nu_prompt_sigma_f);

      std::vector<std::vector<double>> delayed_nu_sigma_f_by_prec;
      if (not H5ReadDataset2D<double>(
            file.Id(), path + "delayed-nu-fission", delayed_nu_sigma_f_by_prec))
      {
        std::vector<double> delayed_nu_sigma_f_1d;
        OpenSnLogicalErrorIf(not H5ReadDataset1D<double>(
                               file.Id(), path + "delayed-nu-fission", delayed_nu_sigma_f_1d),
                             "Failed to read dataset \"" + path + "delayed-nu-fission\".");
        OpenSnLogicalErrorIf(
          mgxs.num_precursors_ != 1 or delayed_nu_sigma_f_1d.size() != mgxs.num_groups_,
          "Dataset \"" + path + "delayed-nu-fission\" has incorrect dimensions.");
        delayed_nu_sigma_f_by_prec.push_back(std::move(delayed_nu_sigma_f_1d));
      }
      std::vector<double> nu_delayed_sigma_f(mgxs.num_groups_, 0.0);
      for (const auto& delayed_family : delayed_nu_sigma_f_by_prec)
        for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
          nu_delayed_sigma_f[g] += delayed_family[g];

      if (nu_prompt_sigma_f.empty())
      {
        nu_prompt_sigma_f = mgxs.nu_sigma_f_;
        for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
          nu_prompt_sigma_f[g] -= nu_delayed_sigma_f[g];
      }

      OpenSnLogicalErrorIf(not IsNonNegative(nu_prompt_sigma_f),
                           "OpenMC delayed fission data imply a negative prompt nu-sigma-f.");
      OpenSnLogicalErrorIf(not IsNonNegative(nu_delayed_sigma_f),
                           "OpenMC delayed nu-sigma-f must be non-negative.");

      std::vector<double> precursor_decay_constants;
      H5ReadDataset1D<double>(file.Id(), path + "decay-rate", precursor_decay_constants);
      OpenSnLogicalErrorIf(precursor_decay_constants.size() != mgxs.num_precursors_,
                           "OpenMC decay-rate data has incorrect size.");

      std::vector<std::vector<double>> chi_delayed_by_prec;
      if (not H5ReadDataset2D<double>(file.Id(), path + "chi-delayed", chi_delayed_by_prec))
      {
        std::vector<double> chi_delayed_1d;
        OpenSnLogicalErrorIf(
          not H5ReadDataset1D<double>(file.Id(), path + "chi-delayed", chi_delayed_1d),
          "Failed to read dataset \"" + path + "chi-delayed\".");
        OpenSnLogicalErrorIf(mgxs.num_precursors_ != 1 or chi_delayed_1d.size() != mgxs.num_groups_,
                             "Dataset \"" + path + "chi-delayed\" has incorrect dimensions.");
        chi_delayed_by_prec.push_back(std::move(chi_delayed_1d));
      }

      std::vector<double> chi_prompt = mgxs.chi_;
      if (H5Has(file.Id(), path + "chi-prompt"))
      {
        H5ReadDataset1D<double>(file.Id(), path + "chi-prompt", chi_prompt);
        OpenSnLogicalErrorIf(chi_prompt.size() != mgxs.num_groups_,
                             "OpenMC chi-prompt data has incorrect size.");
      }
      else
      {
        log.Log0Warning() << "OpenMC delayed-neutron data found in \"" << file_name
                          << "\" without chi-prompt. Using chi as the prompt spectrum.";
      }

      mgxs.nu_prompt_sigma_f_ = std::move(nu_prompt_sigma_f);
      mgxs.nu_delayed_sigma_f_ = std::move(nu_delayed_sigma_f);
      mgxs.precursors_.resize(mgxs.num_precursors_);

      double delayed_total = 0.0;
      for (const auto& delayed_family : delayed_nu_sigma_f_by_prec)
        delayed_total += std::accumulate(delayed_family.begin(), delayed_family.end(), 0.0);

      for (unsigned int j = 0; j < mgxs.num_precursors_; ++j)
      {
        auto& precursor = mgxs.precursors_[j];
        precursor.decay_constant = precursor_decay_constants[j];
        precursor.emission_spectrum = chi_delayed_by_prec[j];
        const double family_total = std::accumulate(
          delayed_nu_sigma_f_by_prec[j].begin(), delayed_nu_sigma_f_by_prec[j].end(), 0.0);
        precursor.fractional_yield = delayed_total > 0.0 ? family_total / delayed_total : 0.0;
      }

      mgxs.production_matrix_.assign(mgxs.num_groups_, std::vector<double>(mgxs.num_groups_, 0.0));
      for (unsigned int gp = 0; gp < mgxs.num_groups_; ++gp)
        for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
          mgxs.production_matrix_[g][gp] = chi_prompt[g] * mgxs.nu_prompt_sigma_f_[gp];
    }

    // Production matrix (computed)
    // TODO: This path uses chi * nu_sigma_f (total production). If delayed-neutron data is
    // introduced here in the future, ensure LBS source/fission-production routines remain
    // consistent with prompt-vs-total production to avoid double counting.
    if (mgxs.production_matrix_.empty())
    {
      mgxs.production_matrix_.resize(mgxs.num_groups_, std::vector<double>(mgxs.num_groups_));
      for (unsigned int gp = 0; gp < mgxs.num_groups_; ++gp)
        for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
          mgxs.production_matrix_[g][gp] = mgxs.chi_[g] * mgxs.nu_sigma_f_[gp];

      OpenSnLogicalErrorIf(
        not IsNonNegative(nu),
        "The production matrix implies an invalid negative average fission neutron yield.");

      OpenSnLogicalErrorIf(
        not std::all_of(nu.begin(), nu.end(), [](double x) { return x == 0.0 or x > 1.0; }),
        "Incompatible fission data encountered. The computed nu is not either zero or "
        "greater than one.");

      if (std::any_of(nu.begin(), nu.end(), [](double x) { return x > 8.0; }))
        log.Log0Warning() << "A computed nu of greater than 8.0 was encountered. ";
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
