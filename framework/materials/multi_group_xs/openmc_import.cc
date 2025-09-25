// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/hdf_utils.h"
#include "framework/utils/utils.h"
#include <numeric>
#include <algorithm>

namespace opensn
{

MultiGroupXS
MultiGroupXS::LoadFromOpenMC(const std::string& file_name,
                             const std::string& dataset_name,
                             double temperature)
{
  MultiGroupXS mgxs;

  // Open file
  hid_t file = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (not file)
  {
    std::string err_msg = "Unable to open " + file_name + " or it is not a valid HDF5 file.\n";
    throw std::logic_error(err_msg);
  }

  // Check file type
  std::string filetype;
  H5ReadAttribute<std::string>(file, "filetype", filetype);
  if (filetype != "mgxs")
    throw std::runtime_error(file_name + " is not a valid OpenMC library file");

  log.Log() << "Reading OpenMC cross-section file \"" << file_name << "\"\n";

  // Number of groups
  if (not H5ReadAttribute<size_t>(file, "energy_groups", mgxs.num_groups_))
    throw std::runtime_error("Failure reading \"energy_groups\" from " + file_name);

  // Group structure
  H5ReadDataset1D<double>(file, "/group structure", mgxs.e_bounds_);
  std::reverse(mgxs.e_bounds_.begin(), mgxs.e_bounds_.end());

  // Temperature
  log.Log0() << file_name + " cross-section data evaluated at " << temperature << "K\n";
  mgxs.temperature_ = temperature;

  // Base path
  std::stringstream ss;
  ss << std::fixed << std::setprecision(0) << mgxs.temperature_ << "K";
  std::string path = "/" + dataset_name + "/" + ss.str() + "/";
  if (!H5Has(file, path))
  {
    throw std::runtime_error("Could not find dataset " + dataset_name + "/" + ss.str() + " in " +
                             file_name);
  }
  mgxs.temperature_ = temperature;

  // Scattering order
  if (not H5ReadGroupAttribute<unsigned int>(file, dataset_name, "order", mgxs.scattering_order_))
    throw std::runtime_error("Failure reading \"order\" from " + file_name);

  // Scattering shape
  std::string scatter_shape;
  H5ReadGroupAttribute<std::string>(file, dataset_name, "scatter_shape", scatter_shape);
  if (scatter_shape != "[G][G'][Order]")
    throw std::runtime_error(file_name + " has an unsupported scatter shape");

  // Number of precursors
  mgxs.num_precursors_ = 0;

  // Inverse velocity
  H5ReadDataset1D<double>(file, path + "inverse-velocity", mgxs.inv_velocity_);
  OpenSnLogicalErrorIf(not IsNonNegative(mgxs.inv_velocity_),
                       "Only positive inverse velocity values are permitted.");

  // Total
  H5ReadDataset1D<double>(file, path + "total", mgxs.sigma_t_);
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
  if (H5Has(file, path + "scatter_data/scatter_matrix"))
  {
    mgxs.transfer_matrices_.assign(mgxs.scattering_order_ + 1,
                                   SparseMatrix(mgxs.num_groups_, mgxs.num_groups_));
    std::vector<double> flat_scatter_matrix;
    H5ReadDataset1D<double>(file, path + "scatter_data/scatter_matrix", flat_scatter_matrix);
    std::vector<int> g_min, g_max;
    H5ReadDataset1D<int>(file, path + "scatter_data/g_min", g_min);
    H5ReadDataset1D<int>(file, path + "scatter_data/g_max", g_max);
    int fidx = 0;
    for (size_t gp = 0; gp < mgxs.num_groups_; ++gp)
      for (int g = g_min[gp]; g <= g_max[gp]; ++g)
        for (unsigned int n = 0; n < mgxs.scattering_order_ + 1; ++n, ++fidx)
          mgxs.transfer_matrices_.at(n).Insert(g - 1, gp, flat_scatter_matrix[fidx]);
  }

  if (mgxs.sigma_a_.empty())
    mgxs.ComputeAbsorption();
  mgxs.ComputeDiffusionParameters();

  // Is fissionable?
  H5ReadGroupAttribute<bool>(file, dataset_name, "fissionable", mgxs.is_fissionable_);
  if (mgxs.is_fissionable_)
  {
    // Fission
    H5ReadDataset1D<double>(file, path + "fission", mgxs.sigma_f_);
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
    H5ReadDataset1D<double>(file, path + "nu-fission", mgxs.nu_sigma_f_);
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
    H5ReadDataset1D<double>(file, path + "chi", mgxs.chi_);
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
    for (size_t g = 0; g < mgxs.num_groups_; ++g)
      if (mgxs.sigma_f_[g] > 0.0)
        nu[g] /= mgxs.sigma_f_[g];

    // Production matrix (computed)
    if (mgxs.production_matrix_.empty())
    {
      mgxs.production_matrix_.resize(mgxs.num_groups_, std::vector<double>(mgxs.num_groups_));
      for (size_t gp = 0; gp < mgxs.num_groups_; ++gp)
        for (size_t g = 0; g < mgxs.num_groups_; ++g)
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
    mgxs.sigma_f_.clear();
    mgxs.nu_sigma_f_.clear();
    mgxs.nu_prompt_sigma_f_.clear();
    mgxs.nu_delayed_sigma_f_.clear();
    mgxs.production_matrix_.clear();
    mgxs.precursors_.clear();
  }

  H5Fclose(file);

  return mgxs;
}

} // namespace opensn
