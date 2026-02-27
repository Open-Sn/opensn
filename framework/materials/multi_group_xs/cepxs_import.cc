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
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace opensn
{
namespace
{

class FortranRecordReader
{
public:
  explicit FortranRecordReader(const std::string& filename) : in_(filename, std::ios::binary) {}

  bool IsOpen() const { return in_.is_open(); }

  bool ReadRecord(std::vector<char>& payload)
  {
    std::uint32_t len = 0;
    if (not in_.read(reinterpret_cast<char*>(&len), sizeof(len)))
      return false;

    payload.resize(len);
    if (len > 0 and not in_.read(payload.data(), static_cast<std::streamsize>(len)))
      throw std::runtime_error("Failed reading Fortran record payload.");

    std::uint32_t tail = 0;
    if (not in_.read(reinterpret_cast<char*>(&tail), sizeof(tail)))
      throw std::runtime_error("Failed reading Fortran record trailer.");
    if (tail != len)
      throw std::runtime_error("Fortran record marker mismatch.");

    return true;
  }

private:
  std::ifstream in_;
};

struct ParsedCEPXSData
{
  unsigned int num_groups = 0;
  unsigned int scattering_order = 0;
  bool is_fissionable = false;
  unsigned int num_precursors = 0;

  std::vector<double> e_bounds;
  std::vector<double> sigma_t;
  std::vector<double> charge_deposition;
  std::vector<double> secondary_production;
  std::vector<double> energy_deposition;
  std::vector<SparseMatrix> transfer_matrices;
};

std::vector<std::int32_t>
BytesToInt32(const std::vector<char>& bytes)
{
  OpenSnLogicalErrorIf(bytes.size() % sizeof(std::int32_t) != 0,
                       "Invalid record size for int32 conversion.");
  std::vector<std::int32_t> vals(bytes.size() / sizeof(std::int32_t), 0);
  std::memcpy(vals.data(), bytes.data(), bytes.size());
  return vals;
}

std::vector<double>
BytesToDouble(const std::vector<char>& bytes)
{
  OpenSnLogicalErrorIf(bytes.size() % sizeof(double) != 0,
                       "Invalid record size for double conversion.");
  std::vector<double> vals(bytes.size() / sizeof(double), 0.0);
  std::memcpy(vals.data(), bytes.data(), bytes.size());
  return vals;
}

std::vector<double>
ExtractEnergyBoundsFromAncillary(const std::vector<char>& ancillary, const int n_groups)
{
  OpenSnLogicalErrorIf(n_groups <= 0, "Invalid group count for CEPXS ancillary parsing.");
  OpenSnLogicalErrorIf(ancillary.size() % sizeof(double) != 0,
                       "CEPXS ancillary record is not an integer multiple of 8 bytes.");

  const auto vals = BytesToDouble(ancillary);
  const size_t n_bounds = static_cast<size_t>(n_groups + 1);
  OpenSnLogicalErrorIf(vals.size() < n_bounds,
                       "CEPXS ancillary record is too short to contain group boundaries.");

  const auto is_valid_bounds = [&](const size_t start_idx)
  {
    const double e0 = vals[start_idx];
    const double eN = vals[start_idx + n_bounds - 1];
    if (not std::isfinite(e0) or not std::isfinite(eN) or e0 <= 0.0 or eN <= 0.0 or e0 <= eN)
      return false;

    for (size_t i = 1; i < n_bounds; ++i)
    {
      const double e_prev = vals[start_idx + i - 1];
      const double e_curr = vals[start_idx + i];
      if (not std::isfinite(e_curr) or e_curr <= 0.0 or e_prev <= e_curr)
        return false;
    }
    return true;
  };

  // CEPXS BFP ancillary records commonly place group boundaries near index 96.
  const size_t canonical_start = 96;
  if (canonical_start + n_bounds <= vals.size() and is_valid_bounds(canonical_start))
    return std::vector<double>(vals.begin() + static_cast<std::ptrdiff_t>(canonical_start),
                               vals.begin() +
                                 static_cast<std::ptrdiff_t>(canonical_start + n_bounds));

  // Fallback: search the entire ancillary vector for a strictly decreasing positive window.
  for (size_t start = 0; start + n_bounds <= vals.size(); ++start)
    if (is_valid_bounds(start))
      return std::vector<double>(vals.begin() + static_cast<std::ptrdiff_t>(start),
                                 vals.begin() + static_cast<std::ptrdiff_t>(start + n_bounds));

  throw std::logic_error("Failed to locate CEPXS group boundaries in binary ancillary record.");
}

bool
LooksLikeFortranBinary(const std::string& filename)
{
  std::ifstream in(filename, std::ios::binary);
  if (not in.is_open())
    return false;

  std::uint32_t marker = 0;
  if (not in.read(reinterpret_cast<char*>(&marker), sizeof(marker)))
    return false;

  return marker > 0 and marker < (1u << 20);
}

ParsedCEPXSData
ParseCEPXSText(const std::string& filename, int material_id)
{
  std::ifstream fin(filename);
  OpenSnLogicalErrorIf(not fin.is_open(), "Unable to open CEPXS file \"" + filename + "\".");
  log.Log() << "Reading CEPXS cross-section file \"" << filename << "\"\n";

  ParsedCEPXSData xs;

  int n_materials = 0;
  std::array<int, 3> n_groups_particle{0, 0, 0}; // gamma, electron, positron
  fin >> n_groups_particle[0] >> n_groups_particle[1] >> n_groups_particle[2] >> n_materials;
  OpenSnLogicalErrorIf(not fin.good(),
                       "Failed to parse CEPXS header in file \"" + filename + "\".");

  OpenSnLogicalErrorIf(n_materials != 1,
                       "CEPXS text reader currently supports exactly one material per file.");
  OpenSnLogicalErrorIf(material_id != 0,
                       "CEPXS text reader only supports material_id=0 for single-material files.");
  OpenSnLogicalErrorIf(std::any_of(n_groups_particle.begin(),
                                   n_groups_particle.end(),
                                   [](int n) { return n < 0; }),
                       "CEPXS group counts must be non-negative.");

  const int num_groups = n_groups_particle[0] + n_groups_particle[1] + n_groups_particle[2];
  OpenSnLogicalErrorIf(num_groups <= 0, "CEPXS file has zero total groups.");
  xs.num_groups = static_cast<unsigned int>(num_groups);

  int L_max = 0;
  fin >> L_max;
  OpenSnLogicalErrorIf(not fin.good(), "Failed parsing CEPXS scattering order.");
  OpenSnLogicalErrorIf(L_max < 0, "CEPXS scattering order must be non-negative.");
  xs.scattering_order = static_cast<unsigned int>(L_max);

  xs.e_bounds.resize(xs.num_groups + 1, 0.0);
  for (unsigned int b = 0; b <= xs.num_groups; ++b)
    xs.e_bounds[b] = static_cast<double>(xs.num_groups - b);

  xs.sigma_t.assign(xs.num_groups, 0.0);
  xs.charge_deposition.clear();
  xs.secondary_production.clear();
  xs.energy_deposition.assign(xs.num_groups, 0.0);
  xs.transfer_matrices.assign(xs.scattering_order + 1, SparseMatrix(xs.num_groups, xs.num_groups));

  auto read_particle_block = [&](std::vector<double>& destination)
  {
    unsigned int offset = 0;
    for (const int n_groups_for_particle : n_groups_particle)
    {
      for (int g = 0; g < n_groups_for_particle; ++g)
      {
        fin >> destination.at(offset + static_cast<unsigned int>(g));
        OpenSnLogicalErrorIf(not fin.good(),
                             "Unexpected end-of-file reading CEPXS block in \"" + filename +
                               "\".");
      }
      offset += static_cast<unsigned int>(n_groups_for_particle);
    }
  };

  read_particle_block(xs.sigma_t);
  read_particle_block(xs.energy_deposition);

  const auto is_finite_vec = [](const std::vector<double>& vec)
  {
    return std::all_of(vec.begin(), vec.end(), [](const double v) { return std::isfinite(v); });
  };

  OpenSnLogicalErrorIf(not IsNonNegative(xs.sigma_t),
                       "CEPXS total cross section contains negative values.");
  OpenSnLogicalErrorIf(not is_finite_vec(xs.sigma_t),
                       "CEPXS total cross section contains non-finite values.");
  OpenSnLogicalErrorIf(not is_finite_vec(xs.energy_deposition),
                       "CEPXS energy deposition contains non-finite values.");

  for (unsigned int ell = 0; ell <= xs.scattering_order; ++ell)
  {
    auto& Sm = xs.transfer_matrices[ell];
    for (unsigned int g = 0; g < xs.num_groups; ++g)
      for (unsigned int gp = 0; gp < xs.num_groups; ++gp)
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

  return xs;
}

ParsedCEPXSData
ParseCEPXSBFPBinary(const std::string& filename, int material_id)
{
  FortranRecordReader rdr(filename);
  OpenSnLogicalErrorIf(not rdr.IsOpen(),
                       "Unable to open CEPXS binary file \"" + filename + "\".");
  log.Log() << "Reading CEPXS-BFP binary cross-section file \"" << filename << "\"\n";

  ParsedCEPXSData xs;

  std::vector<char> rec;
  OpenSnLogicalErrorIf(not rdr.ReadRecord(rec), "Failed reading CEPXS binary title record.");

  OpenSnLogicalErrorIf(not rdr.ReadRecord(rec), "Failed reading CEPXS binary metadata record.");
  const auto meta = BytesToInt32(rec);
  OpenSnLogicalErrorIf(meta.size() < 8, "CEPXS binary metadata record is too short.");

  const int n_groups = meta[0];
  const int n_materials = meta[1];
  const int n_entries = meta[2];
  const int total_xs_row = meta[3] - 1;      // Convert to 0-based indexing.
  const int self_scatter_row = meta[4] - 1;  // Convert to 0-based indexing.
  const int n_moments = meta[5];
  const int n_tables_from_header = meta[7];

  OpenSnLogicalErrorIf(n_materials <= 0, "CEPXS binary has invalid number of materials.");
  OpenSnLogicalErrorIf(n_groups <= 0, "CEPXS binary has invalid number of groups.");
  OpenSnLogicalErrorIf(n_entries <= 8, "CEPXS binary has invalid number of entries.");
  OpenSnLogicalErrorIf(total_xs_row < 0 || total_xs_row >= n_entries,
                       "CEPXS binary has invalid total-xs row index.");
  OpenSnLogicalErrorIf(self_scatter_row < 0 || self_scatter_row >= n_entries,
                       "CEPXS binary has invalid self-scatter row index.");
  OpenSnLogicalErrorIf(n_moments <= 0, "CEPXS binary has invalid number of moments.");
  OpenSnLogicalErrorIf(material_id < 0 || material_id >= n_materials,
                       "CEPXS binary material_id out of range.");

  OpenSnLogicalErrorIf(not rdr.ReadRecord(rec), "Failed reading CEPXS binary ancillary record.");

  xs.num_groups = static_cast<unsigned int>(n_groups);
  xs.e_bounds = ExtractEnergyBoundsFromAncillary(rec, n_groups);

  xs.sigma_t.assign(xs.num_groups, 0.0);
  xs.charge_deposition.assign(xs.num_groups, 0.0);
  xs.secondary_production.assign(xs.num_groups, 0.0);
  xs.energy_deposition.assign(xs.num_groups, 0.0);

  std::vector<std::vector<double>> moment_tables;
  while (rdr.ReadRecord(rec))
  {
    OpenSnLogicalErrorIf(rec.size() != static_cast<size_t>(n_groups * n_entries * sizeof(double)),
                         "Unexpected CEPXS binary moment-record size.");
    moment_tables.push_back(BytesToDouble(rec));
  }

  OpenSnLogicalErrorIf(moment_tables.empty(), "CEPXS binary contains no moment records.");
  if (n_tables_from_header > 0)
    OpenSnLogicalErrorIf(static_cast<int>(moment_tables.size()) != n_tables_from_header,
                         "CEPXS binary table count mismatch with header.");

  OpenSnLogicalErrorIf(static_cast<int>(moment_tables.size()) != n_materials * n_moments,
                       "CEPXS binary table count does not match materials*moments.");

  xs.scattering_order = static_cast<unsigned int>(n_moments - 1);
  xs.transfer_matrices.assign(xs.scattering_order + 1, SparseMatrix(xs.num_groups, xs.num_groups));
  constexpr int charge_deposition_row = 0; // 1-based row 1
  constexpr int secondary_production_row = 1; // 1-based row 2
  constexpr int energy_deposition_row = 2; // 1-based row 3
  const int first_transfer_row = std::min(self_scatter_row, total_xs_row + 1);

  for (int mom = 0; mom < n_moments; ++mom)
  {
    const int table_idx = material_id * n_moments + mom;
    OpenSnLogicalErrorIf(table_idx < 0 ||
                           table_idx >= static_cast<int>(moment_tables.size()),
                         "Computed CEPXS binary table index out of range.");

    const auto& table = moment_tables[table_idx];
    auto& Sm = xs.transfer_matrices[static_cast<size_t>(mom)];

    for (int g_to = 0; g_to < n_groups; ++g_to)
      for (int row = 0; row < n_entries; ++row)
      {
        const double value = table[static_cast<size_t>(row + n_entries * g_to)];

        if (mom == 0)
        {
          if (row == charge_deposition_row)
            xs.charge_deposition[g_to] = value;
          else if (row == secondary_production_row)
            xs.secondary_production[g_to] = value;
          else if (row == energy_deposition_row)
            xs.energy_deposition[g_to] = value;
          else if (row == total_xs_row)
            xs.sigma_t[g_to] = value;
        }

        if (row < first_transfer_row || value == 0.0)
          continue;

        int g_from = -1;
        if (row < self_scatter_row)
          g_from = self_scatter_row - row + g_to;
        else if (row == self_scatter_row)
          g_from = g_to;
        else
          g_from = g_to - (row - self_scatter_row);

        if (g_from < 0 or g_from >= n_groups)
          continue;

        Sm.Insert(g_to, g_from, value);
      }
  }

  const auto is_finite_vec = [](const std::vector<double>& vec)
  {
    return std::all_of(vec.begin(), vec.end(), [](const double v) { return std::isfinite(v); });
  };

  OpenSnLogicalErrorIf(not IsNonNegative(xs.sigma_t),
                       "CEPXS binary total cross section contains negative values.");
  OpenSnLogicalErrorIf(not is_finite_vec(xs.sigma_t),
                       "CEPXS binary total cross section contains non-finite values.");
  OpenSnLogicalErrorIf(not is_finite_vec(xs.energy_deposition),
                       "CEPXS binary energy deposition contains non-finite values.");

  return xs;
}

} // namespace

MultiGroupXS
MultiGroupXS::LoadFromCEPXS(const std::string& filename, int material_id)
{
  MultiGroupXS mgxs;
  const auto parsed = LooksLikeFortranBinary(filename) ? ParseCEPXSBFPBinary(filename, material_id)
                                                        : ParseCEPXSText(filename, material_id);

  mgxs.num_groups_ = parsed.num_groups;
  mgxs.scattering_order_ = parsed.scattering_order;
  mgxs.is_fissionable_ = parsed.is_fissionable;
  mgxs.num_precursors_ = parsed.num_precursors;

  mgxs.e_bounds_ = parsed.e_bounds;
  mgxs.sigma_t_ = parsed.sigma_t;
  // Derive absorption from total and transfer matrices for consistency with
  // the rest of OpenSn's XS handling.
  mgxs.sigma_a_.clear();
  mgxs.energy_deposition_ = parsed.energy_deposition;
  mgxs.transfer_matrices_ = parsed.transfer_matrices;
  if (not parsed.charge_deposition.empty())
    mgxs.custom_xs_["cepxs_charge_deposition"] = parsed.charge_deposition;
  if (not parsed.secondary_production.empty())
    mgxs.custom_xs_["cepxs_secondary_production"] = parsed.secondary_production;

  mgxs.ComputeAbsorption();
  mgxs.ComputeDiffusionParameters();
  return mgxs;
}

} // namespace opensn
