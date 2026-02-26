// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"

#include <fstream>
#include <array>
#include <algorithm>
#include <cmath>

namespace opensn
{

MultiGroupXS
MultiGroupXS::LoadFromCEPXS(const std::string& filename)
{
  std::ifstream fin(filename);
  OpenSnLogicalErrorIf(not fin.is_open(), "Unable to open CEPXS file \"" + filename + "\".");
  log.Log() << "Reading CEPXS cross-section file \"" << filename << "\"\n";

  MultiGroupXS mgxs;

  int n_materials = 0;
  std::array<int, 3> n_groups_particle{0, 0, 0}; // gamma, electron, positron
  fin >> n_groups_particle[0] >> n_groups_particle[1] >> n_groups_particle[2] >> n_materials;
  OpenSnLogicalErrorIf(not fin.good(),
                       "Failed to parse CEPXS header in file \"" + filename + "\".");

  OpenSnLogicalErrorIf(n_materials != 1,
                       "CEPXS reader currently supports exactly one material per file.");
  OpenSnLogicalErrorIf(std::any_of(n_groups_particle.begin(),
                                   n_groups_particle.end(),
                                   [](int n) { return n < 0; }),
                       "CEPXS group counts must be non-negative.");

  const int num_groups = n_groups_particle[0] + n_groups_particle[1] + n_groups_particle[2];
  OpenSnLogicalErrorIf(num_groups <= 0, "CEPXS file has zero total groups.");
  mgxs.num_groups_ = static_cast<unsigned int>(num_groups);

  int L_max = 0;
  fin >> L_max;
  OpenSnLogicalErrorIf(not fin.good(), "Failed parsing CEPXS scattering order.");
  OpenSnLogicalErrorIf(L_max < 0, "CEPXS scattering order must be non-negative.");
  mgxs.scattering_order_ = static_cast<unsigned int>(L_max);

  mgxs.is_fissionable_ = false;
  mgxs.num_precursors_ = 0;

  // CEPXS does not provide native group boundaries in this form. We maintain a monotonic
  // descending pseudo-boundary set for introspection consistency.
  mgxs.e_bounds_.resize(mgxs.num_groups_ + 1, 0.0);
  for (unsigned int b = 0; b <= mgxs.num_groups_; ++b)
    mgxs.e_bounds_[b] = static_cast<double>(mgxs.num_groups_ - b);

  mgxs.sigma_t_.assign(mgxs.num_groups_, 0.0);
  mgxs.energy_deposition_.assign(mgxs.num_groups_, 0.0);
  mgxs.transfer_matrices_.assign(mgxs.scattering_order_ + 1,
                                 SparseMatrix(mgxs.num_groups_, mgxs.num_groups_));

  auto read_particle_block = [&](std::vector<double>& destination)
  {
    unsigned int offset = 0;
    for (const int n_groups_for_particle : n_groups_particle)
    {
      for (int g = 0; g < n_groups_for_particle; ++g)
      {
        fin >> destination.at(offset + static_cast<unsigned int>(g));
        OpenSnLogicalErrorIf(not fin.good(),
                             "Unexpected end-of-file reading CEPXS block in \"" + filename + "\".");
      }
      offset += static_cast<unsigned int>(n_groups_for_particle);
    }
  };

  read_particle_block(mgxs.sigma_t_);
  read_particle_block(mgxs.energy_deposition_);

  const auto is_finite_vec = [](const std::vector<double>& vec)
  {
    return std::all_of(vec.begin(), vec.end(), [](const double v) { return std::isfinite(v); });
  };

  OpenSnLogicalErrorIf(not IsNonNegative(mgxs.sigma_t_),
                       "CEPXS total cross section contains negative values.");
  OpenSnLogicalErrorIf(not is_finite_vec(mgxs.sigma_t_),
                       "CEPXS total cross section contains non-finite values.");
  OpenSnLogicalErrorIf(not is_finite_vec(mgxs.energy_deposition_),
                       "CEPXS energy deposition contains non-finite values.");

  for (unsigned int ell = 0; ell <= mgxs.scattering_order_; ++ell)
  {
    auto& Sm = mgxs.transfer_matrices_[ell];
    for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
      for (unsigned int gp = 0; gp < mgxs.num_groups_; ++gp)
      {
        double val = 0.0;
        fin >> val;
        OpenSnLogicalErrorIf(not fin.good(),
                             "Unexpected end-of-file reading CEPXS transfer matrices in \"" +
                               filename + "\".");
        OpenSnLogicalErrorIf(not std::isfinite(val),
                             "CEPXS transfer matrix contains non-finite values.");
        if (val != 0.0)
          Sm.Insert(g, gp, val);
      }
  }

  mgxs.ComputeAbsorption();
  mgxs.ComputeDiffusionParameters();

  return mgxs;
}

} // namespace opensn
