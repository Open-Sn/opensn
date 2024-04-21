// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/logging/log.h"
#include "framework/utils/hdf_utils.h"
#include <numeric>
#include <algorithm>

namespace opensn
{

void
MultiGroupXS::Initialize(const std::string& file_name, double temperature)
{
  Reset();

  // Disable internal HDF error reporting
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  // Open file
  try
  {
    H5::H5File::isHdf5(file_name.c_str());
  }
  catch (H5::FileIException(&err))
  {
    std::string err_msg = file_name + " is not an HDF5 file.\n";
    throw std::logic_error(err_msg);
  }
  H5::H5File file(file_name.c_str(), H5F_ACC_RDONLY);

  // Check file type
  std::string filetype;
  H5ReadAttribute<std::string>(file, "filetype", filetype);
  if (filetype != "mgxs")
    throw std::runtime_error(file_name + " is not a valid OpenMC library file");

  log.Log() << "Reading OpenMC cross section file \"" << file_name << "\"\n";

  // Number of groups
  if (not H5ReadAttribute<size_t>(file, "energy_groups", num_groups_))
    throw std::runtime_error("Failure reading \"energy_groups\" from " + file_name);

  // Group structure
  e_bounds_ = H5ReadDataset1D<double>(file, "/group structure");
  std::reverse(e_bounds_.begin(), e_bounds_.end());

  // Temperature
  log.Log0() << file_name + " cross-section data evaluated at " << temperature << "K\n";
  temperature_ = temperature;

  // Base path
  std::stringstream ss;
  ss << std::fixed << std::setprecision(0) << temperature_ << "K";
  std::string path = "/set1/" + ss.str() + "/";
  if (!H5Has(file, path))
    throw std::runtime_error("Could not find data for " + ss.str() + " in " + file_name);
  temperature_ = temperature;

  // Scattering order
  if (not H5ReadGroupAttribute<size_t>(file, "set1", "order", scattering_order_))
    throw std::runtime_error("Failure reading \"order\" from " + file_name);

  // Number of precursors
  num_precursors_ = 0;

  // Velocity
  inv_velocity_ = H5ReadDataset1D<double>(file, path + "inverse-velocity");
  OpenSnLogicalErrorIf(not IsNonNegative(inv_velocity_),
                       "Only positive inverse velocity values are permitted.");

  // Total
  sigma_t_ = H5ReadDataset1D<double>(file, path + "total");
  OpenSnLogicalErrorIf(sigma_t_.empty(),
                       "\"total\" data block not found in " + file_name + ".");
  OpenSnLogicalErrorIf(not IsNonNegative(sigma_t_),
                       "Only non-negative total cross section values are permitted.");

  // Absorption
  sigma_a_ = H5ReadDataset1D<double>(file, path + "absorption");
  OpenSnLogicalErrorIf(not IsNonNegative(sigma_a_),
                       "Only non-negative absorption cross section values are permitted.");

  // Transfer
  if (H5Has(file, path + "scatter_data/scatter_matrix"))
  {
    transfer_matrices_.assign(scattering_order_ + 1, SparseMatrix(num_groups_, num_groups_));
    auto flat_scatter_matrix = H5ReadDataset1D<double>(file, path + "scatter_data/scatter_matrix");
    auto g_min = H5ReadDataset1D<int>(file, path + "scatter_data/g_min");
    auto g_max = H5ReadDataset1D<int>(file, path + "scatter_data/g_max");
    int fidx = 0;
    for (int g = 0; g < num_groups_; ++g)
      for (int gp = g_min[g]; gp <= g_max[g]; ++gp)
        for (int n = 0; n < scattering_order_ + 1; ++n, ++fidx)
          transfer_matrices_.at(n).Insert(g, gp - 1, flat_scatter_matrix[fidx]);
  }

  if (sigma_a_.empty())
    ComputeAbsorption();
  ComputeDiffusionParameters();

  // Is fissionable?
  H5ReadGroupAttribute<bool>(file, "set1", "fissionable", is_fissionable_);
  if (is_fissionable_)
  {
    // Fission
    sigma_f_ = H5ReadDataset1D<double>(file, path + "fission");
    OpenSnLogicalErrorIf(sigma_f_.empty(),
                         "\"fission\" data block not found in " + file_name + ".");
    OpenSnLogicalErrorIf(not IsNonNegative(sigma_f_),
                         "Only non-negative fission cross section values are permitted.");
    if (not HasNonZero(sigma_f_))
    {
      log.Log0Warning() << "The fission cross section specified in "
                        << "\"" << file_name << "\" is uniformly zero... Clearing it.";
      sigma_f_.clear();
    }

    // Nu-Fission
    nu_sigma_f_ = H5ReadDataset1D<double>(file, path + "nu-fission");
    OpenSnLogicalErrorIf(nu_sigma_f_.empty(),
                         "\"nu-fission\" data block not found in " + file_name + ".");
    OpenSnLogicalErrorIf(
      not IsNonNegative(nu_sigma_f_),
      "Only non-negative total fission multiplication cross section values are permitted.");
    if (not HasNonZero(nu_sigma_f_))
    {
      log.Log0Warning() << "The production cross-section specified in "
                        << "\"" << file_name << "\" is uniformly zero... Clearing it.";
      nu_sigma_f_.clear();
    }

    // Chi
    auto chi = H5ReadDataset1D<double>(file, path + "chi");
    OpenSnLogicalErrorIf(chi.empty(),
                         "\"chi\" data block not found in " + file_name + ".");
    if (not chi.empty())
    {
      OpenSnLogicalErrorIf(not HasNonZero(chi),
                           "Steady-state fission spectrum must have at least one non-zero value.");
      OpenSnLogicalErrorIf(not IsNonNegative(chi),
                           "Steady-state fission spectrum must be non-negative.");
      // Normalizing
      const auto sum = std::accumulate(chi.begin(), chi.end(), 0.0);
      std::transform(chi.begin(), chi.end(), chi.begin(), [sum](double& x) { return x / sum; });
    }

    // Nu (computed)
    auto nu = nu_sigma_f_;
    for (size_t g = 0; g < num_groups_; ++g)
      if (sigma_f_[g] > 0.0)
        nu[g] /= sigma_f_[g];

    // Production matrix (computed)
    if (production_matrix_.empty())
    {
      production_matrix_.resize(num_groups_, std::vector<double>(num_groups_));
      for (size_t g = 0; g < num_groups_; ++g)
        for (size_t gp = 0; gp < num_groups_; ++gp)
          production_matrix_[g][gp] = chi[gp] * nu_sigma_f_[g];

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
    sigma_f_.clear();
    nu_sigma_f_.clear();
    nu_prompt_sigma_f_.clear();
    nu_delayed_sigma_f_.clear();
    production_matrix_.clear();
    precursors_.clear();
  }
}

} // namespace opensn
