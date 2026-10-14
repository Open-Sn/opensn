// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

/*
 * CEPXS uses the BXSLIB binary cross-section file format. A description of the
 * BXSLIB format can be found in:
 *
 *   Shapiro, A., and Huria, H.
 *   "Standard Interface File Format"
 *   University of Cincinnati Nuclear Engineering Program, 1992
 *   https://www.osti.gov/servlets/purl/10115350
 *
 * CEPXS format blocks for various materials can be found in:
 *
 *   McConn, R. J., Gesh, C. J., Pagh, R. T., et al.
 *   "Compendium of Material Composition Data for Radiation Transport Modeling"
 *   Pacific Northwest National Laboratory, PNNL-15870 Rev. 1, 2006
 *   https://www.pnnl.gov/main/publications/external/technical_reports/
 *   pnnl-15870rev1.pdf
 */

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include <fstream>
#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace opensn
{
namespace
{

constexpr std::size_t MAX_CEPXS_RECORD_BYTES = std::size_t{256} * 1024U * 1024U;

class FortranRecordReader
{
public:
  explicit FortranRecordReader(const std::string& filename) : in_(filename, std::ios::binary) {}

  bool IsOpen() const { return in_.is_open(); }

  bool ReadRecord(std::vector<char>& payload,
                  const std::size_t max_payload_bytes = MAX_CEPXS_RECORD_BYTES)
  {
    std::uint32_t len = 0;
    if (not ReadU32(in_, len, true))
      return false;

    if (static_cast<std::size_t>(len) > max_payload_bytes)
      throw std::runtime_error("CEPXS Fortran record exceeds the permitted size.");

    payload.resize(len);
    if (len > 0 and not in_.read(payload.data(), static_cast<std::streamsize>(len)))
      throw std::runtime_error("Failed reading Fortran record payload.");

    std::uint32_t tail = 0;
    ReadU32(in_, tail, false);
    if (tail != len)
      throw std::runtime_error("Fortran record marker mismatch.");

    return true;
  }

private:
  static bool ReadU32(std::istream& in, std::uint32_t& value, const bool allow_clean_eof)
  {
    std::array<char, sizeof(std::uint32_t)> bytes{};
    if (not in.read(bytes.data(), static_cast<std::streamsize>(bytes.size())))
    {
      if (allow_clean_eof and in.eof() and in.gcount() == 0)
        return false;
      throw std::runtime_error(allow_clean_eof ? "Failed reading Fortran record header."
                                               : "Failed reading Fortran record trailer.");
    }
    std::memcpy(&value, bytes.data(), sizeof(value));
    return true;
  }

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
  std::vector<double> stopping_power;
  std::vector<SparseMatrix> transfer_matrices;
};

using GroupRange = std::pair<std::size_t, std::size_t>;

std::vector<GroupRange>
FindParticleGroupRanges(const std::vector<double>& bounds)
{
  OpenSnLogicalErrorIf(bounds.size() < 2, "CEPXS energy group structure is incomplete.");

  std::vector<GroupRange> ranges;
  std::size_t begin = 0;
  const std::size_t num_groups = bounds.size() - 1;
  for (std::size_t g = 1; g < num_groups; ++g)
    if (bounds[g] <= bounds[g + 1])
    {
      ranges.emplace_back(begin, g);
      begin = g;
    }
  ranges.emplace_back(begin, num_groups);
  return ranges;
}

ParticleType
ParseParticleType(std::string name)
{
  std::transform(name.begin(),
                 name.end(),
                 name.begin(),
                 [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (name == "photon")
    return ParticleType::PHOTON;
  if (name == "electron")
    return ParticleType::ELECTRON;
  if (name == "positron")
    return ParticleType::POSITRON;
  throw std::invalid_argument("Invalid CEPXS particle type \"" + name +
                              "\". Expected photon, electron, or positron.");
}

std::vector<ParticleType>
ResolveParticleTypes(const std::vector<GroupRange>& ranges,
                     const std::vector<double>& stopping_power,
                     const std::vector<std::string>& particle_order)
{
  std::vector<ParticleType> block_types(ranges.size(), ParticleType::UNKNOWN);
  const bool has_stopping_power = not stopping_power.empty();
  const auto has_nonzero_stopping = [&stopping_power](const GroupRange& range)
  {
    for (std::size_t g = range.first; g < range.second; ++g)
      if (stopping_power[g] > MultiGroupXS::STOPPING_POWER_TOLERANCE)
        return true;
    return false;
  };

  if (not particle_order.empty())
  {
    OpenSnInvalidArgumentIf(particle_order.size() != ranges.size(),
                            "CEPXS particle_order must contain one entry for each of the " +
                              std::to_string(ranges.size()) +
                              " particle blocks found in the BXSLIB energy structure.");
    std::transform(
      particle_order.begin(), particle_order.end(), block_types.begin(), ParseParticleType);
    for (std::size_t i = 0; i < block_types.size(); ++i)
      for (std::size_t j = i + 1; j < block_types.size(); ++j)
        OpenSnInvalidArgumentIf(block_types[i] == block_types[j],
                                "CEPXS particle_order supports at most one block of each particle "
                                "type.");
  }
  else if (has_stopping_power and ranges.size() == 1)
  {
    if (has_nonzero_stopping(ranges.front()))
      block_types.front() = ParticleType::ELECTRON;
  }

  OpenSnInvalidArgumentIf(has_stopping_power and ranges.size() > 1 and particle_order.empty(),
                          "A multi-block CEPXS CSDA library requires particle_order.");

  std::vector<ParticleType> particle_types(ranges.back().second, ParticleType::UNKNOWN);
  for (std::size_t b = 0; b < ranges.size(); ++b)
  {
    const auto [begin, end] = ranges[b];
    if (has_stopping_power and block_types[b] != ParticleType::UNKNOWN)
    {
      OpenSnInvalidArgumentIf(block_types[b] == ParticleType::PHOTON and
                                has_nonzero_stopping(ranges[b]),
                              "CEPXS particle_order labels a block with nonzero stopping power as "
                              "photon.");
    }
    for (std::size_t g = begin; g < end; ++g)
      particle_types[g] = block_types[b];
  }
  return particle_types;
}

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
ExtractEnergyGroupStructure(const std::vector<char>& group_structure_record,
                            const int n_groups,
                            const bool require_energy_bounds)
{
  OpenSnLogicalErrorIf(n_groups <= 0, "Invalid group count for CEPXS group-structure parsing.");
  OpenSnLogicalErrorIf(group_structure_record.size() % sizeof(double) != 0,
                       "CEPXS group-structure record is not an integer multiple of 8 bytes.");

  const auto vals = BytesToDouble(group_structure_record);
  const auto n_bounds = static_cast<size_t>(n_groups) + 1U;
  OpenSnLogicalErrorIf(vals.size() < n_bounds,
                       "CEPXS group-structure record is too short to contain group boundaries.");

  const auto is_valid_bounds = [&](const size_t start_idx)
  {
    const double e0 = vals[start_idx];
    const double eN = vals[start_idx + n_bounds - 1];
    if (not std::isfinite(e0) or not std::isfinite(eN) or e0 <= 0.0 or eN < 0.0 or e0 <= eN)
      return false;

    for (size_t i = 1; i < n_bounds; ++i)
    {
      const double e_prev = vals[start_idx + i - 1];
      const double e_curr = vals[start_idx + i];
      const bool is_terminal_bound = (i + 1 == n_bounds);
      if (not std::isfinite(e_curr) or (is_terminal_bound ? e_curr < 0.0 : e_curr <= 0.0) or
          e_prev <= e_curr)
        return false;
    }
    return true;
  };

  const auto is_finite_positive_window = [&](const size_t start_idx)
  {
    bool any_change = false;
    for (size_t i = 0; i < n_bounds; ++i)
    {
      const double e = vals[start_idx + i];
      if (not std::isfinite(e) or e <= 0.0)
        return false;
      if (i > 0 and vals[start_idx + i - 1] != e)
        any_change = true;
    }
    return any_change;
  };

  // CEPXS BFP group-structure records commonly place group boundaries near index 96.
  const size_t canonical_start = 96;
  if (canonical_start + n_bounds <= vals.size() and is_valid_bounds(canonical_start))
    return {vals.begin() + static_cast<std::ptrdiff_t>(canonical_start),
            vals.begin() + static_cast<std::ptrdiff_t>(canonical_start + n_bounds)};

  // Fallback: search the entire group-structure record for a strictly decreasing positive window.
  for (size_t start = 0; start + n_bounds <= vals.size(); ++start)
    if (is_valid_bounds(start))
      return {vals.begin() + static_cast<std::ptrdiff_t>(start),
              vals.begin() + static_cast<std::ptrdiff_t>(start + n_bounds)};

  // Coupled electron-photon libraries can carry group windows. Use a looser window
  // search before giving up.
  if (canonical_start + n_bounds <= vals.size() and is_finite_positive_window(canonical_start))
    return {vals.begin() + static_cast<std::ptrdiff_t>(canonical_start),
            vals.begin() + static_cast<std::ptrdiff_t>(canonical_start + n_bounds)};

  for (size_t start = 0; start + n_bounds <= vals.size(); ++start)
    if (is_finite_positive_window(start))
      return {vals.begin() + static_cast<std::ptrdiff_t>(start),
              vals.begin() + static_cast<std::ptrdiff_t>(start + n_bounds)};

  OpenSnLogicalErrorIf(require_energy_bounds,
                       "CEPXS energy group structure was not found in the group-structure record. "
                       "CSDA requires physical energy group bounds.");

  // Last resort
  std::vector<double> synthetic(n_bounds, 0.0);
  for (size_t i = 0; i < n_bounds; ++i)
    synthetic[i] = static_cast<double>(n_groups - static_cast<int>(i));
  log.Log()
    << "Warning: CEPXS energy group structure was not found in the group-structure record; using "
       "synthetic bounds.\n";
  return synthetic;
}

bool
LooksLikeFortranBinary(const std::string& filename)
{
  std::ifstream in(filename, std::ios::binary);
  if (not in.is_open())
    return false;

  std::uint32_t marker = 0;
  {
    std::array<char, sizeof(std::uint32_t)> bytes{};
    if (not in.read(bytes.data(), static_cast<std::streamsize>(bytes.size())))
      return false;
    std::memcpy(&marker, bytes.data(), sizeof(marker));
  }

  return marker > 0 and marker < (1U << 20);
}

enum class CEPXSRowFormat
{
  LEGACY = 0,
  CSDA = 1
};

ParsedCEPXSData
ParseCEPXSBFPBinary(const std::string& filename, int material_id, CEPXSRowFormat row_format)
{
  FortranRecordReader rdr(filename);
  OpenSnLogicalErrorIf(not rdr.IsOpen(), "Unable to open CEPXS binary file \"" + filename + "\".");
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
  int total_xs_row = meta[3] - 1;     // Convert to 0-based indexing.
  int self_scatter_row = meta[4] - 1; // Convert to 0-based indexing.
  const int n_moments = meta[5];
  const int n_tables_from_header = meta[7];

  OpenSnLogicalErrorIf(n_materials <= 0, "CEPXS binary has invalid number of materials.");
  OpenSnLogicalErrorIf(n_groups <= 0, "CEPXS binary has invalid number of groups.");
  OpenSnLogicalErrorIf(n_entries <= 8, "CEPXS binary has invalid number of entries.");
  if (row_format == CEPXSRowFormat::LEGACY)
  {
    OpenSnLogicalErrorIf(total_xs_row < 0 or total_xs_row >= n_entries,
                         "CEPXS binary has invalid total-xs row index.");
    OpenSnLogicalErrorIf(self_scatter_row < 0 or self_scatter_row >= n_entries,
                         "CEPXS binary has invalid self-scatter row index.");
  }
  OpenSnLogicalErrorIf(n_moments <= 0, "CEPXS binary has invalid number of moments.");
  OpenSnLogicalErrorIf(material_id < 0 || material_id >= n_materials,
                       "CEPXS binary material_id out of range.");

  const auto num_groups_size = static_cast<std::size_t>(n_groups);
  const auto num_entries_size = static_cast<std::size_t>(n_entries);
  OpenSnLogicalErrorIf(num_groups_size > std::numeric_limits<std::size_t>::max() /
                                           num_entries_size / sizeof(double),
                       "CEPXS binary moment-record size overflows size_t.");
  const auto expected_record_size = num_groups_size * num_entries_size * sizeof(double);
  OpenSnLogicalErrorIf(expected_record_size > MAX_CEPXS_RECORD_BYTES,
                       "CEPXS binary moment record exceeds the supported size limit.");

  const auto num_materials_size = static_cast<std::size_t>(n_materials);
  const auto num_moments_size = static_cast<std::size_t>(n_moments);
  OpenSnLogicalErrorIf(num_materials_size >
                         std::numeric_limits<std::size_t>::max() / num_moments_size,
                       "CEPXS binary table count overflows size_t.");
  const auto expected_num_tables = num_materials_size * num_moments_size;

  OpenSnLogicalErrorIf(not rdr.ReadRecord(rec),
                       "Failed reading CEPXS binary group-structure record.");

  xs.num_groups = static_cast<unsigned int>(n_groups);
  xs.e_bounds = ExtractEnergyGroupStructure(rec, n_groups, row_format == CEPXSRowFormat::CSDA);
  const auto particle_ranges = FindParticleGroupRanges(xs.e_bounds);

  xs.sigma_t.assign(xs.num_groups, 0.0);
  xs.charge_deposition.assign(xs.num_groups, 0.0);
  xs.secondary_production.assign(xs.num_groups, 0.0);
  xs.energy_deposition.assign(xs.num_groups, 0.0);
  // Legacy libraries have no stopping-power row; leave it empty rather than all zeros.
  if (row_format == CEPXSRowFormat::CSDA)
    xs.stopping_power.assign(xs.num_groups, 0.0);

  std::vector<std::vector<double>> moment_tables;
  for (std::size_t table = 0; table < expected_num_tables; ++table)
  {
    OpenSnLogicalErrorIf(not rdr.ReadRecord(rec, expected_record_size),
                         "CEPXS binary contains fewer moment records than declared.");
    OpenSnLogicalErrorIf(rec.size() != expected_record_size,
                         "Unexpected CEPXS binary moment-record size.");
    moment_tables.push_back(BytesToDouble(rec));
  }
  OpenSnLogicalErrorIf(rdr.ReadRecord(rec),
                       "CEPXS binary contains more moment records than declared.");

  OpenSnLogicalErrorIf(moment_tables.empty(), "CEPXS binary contains no moment records.");
  if (n_tables_from_header > 0)
    OpenSnLogicalErrorIf(moment_tables.size() != static_cast<size_t>(n_tables_from_header),
                         "CEPXS binary table count mismatch with header.");

  OpenSnLogicalErrorIf(moment_tables.size() != expected_num_tables,
                       "CEPXS binary table count does not match materials*moments.");

  xs.scattering_order = static_cast<unsigned int>(n_moments - 1);
  xs.transfer_matrices.assign(xs.scattering_order + 1, SparseMatrix(xs.num_groups, xs.num_groups));
  int charge_deposition_row = 0;    // 1-based row 1
  int secondary_production_row = 1; // 1-based row 2
  int energy_deposition_row = 2;    // 1-based row 3
  int stopping_power_row = -1;

  if (row_format == CEPXSRowFormat::CSDA)
  {
    secondary_production_row = 0; // row 1
    charge_deposition_row = 1;    // row 2
    energy_deposition_row = 2;    // row 3
    stopping_power_row = 4;       // row 5
  }

  OpenSnLogicalErrorIf(total_xs_row < 0 or total_xs_row >= n_entries,
                       "CEPXS binary has invalid total-xs row index.");
  OpenSnLogicalErrorIf(self_scatter_row < 0 or self_scatter_row >= n_entries,
                       "CEPXS binary has invalid self-scatter row index.");
  if (stopping_power_row >= 0)
  {
    OpenSnLogicalErrorIf(stopping_power_row >= n_entries,
                         "CEPXS binary has invalid stopping-power row index.");
  }

  const int first_transfer_row = std::min(self_scatter_row, total_xs_row + 1);

  for (int mom = 0; mom < n_moments; ++mom)
  {
    const int table_idx = material_id * n_moments + mom;
    OpenSnLogicalErrorIf(table_idx < 0 || static_cast<size_t>(table_idx) >= moment_tables.size(),
                         "Computed CEPXS binary table index out of range.");

    const auto& table = moment_tables[table_idx];
    auto& Sm = xs.transfer_matrices[static_cast<size_t>(mom)];

    for (int g_to = 0; g_to < n_groups; ++g_to)
      for (int row = 0; row < n_entries; ++row)
      {
        const auto table_index =
          static_cast<size_t>(row) + static_cast<size_t>(n_entries) * static_cast<size_t>(g_to);
        const double value = table[table_index];

        if (mom == 0)
        {
          if (row == charge_deposition_row)
            xs.charge_deposition[g_to] = value;
          else if (row == secondary_production_row)
            xs.secondary_production[g_to] = value;
          else if (row == energy_deposition_row)
            xs.energy_deposition[g_to] = value;
          else if (row == stopping_power_row)
            xs.stopping_power[g_to] = value;
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
  { return std::all_of(vec.begin(), vec.end(), [](const double v) { return std::isfinite(v); }); };

  OpenSnLogicalErrorIf(not IsNonNegative(xs.sigma_t),
                       "CEPXS binary total cross section contains negative values.");
  OpenSnLogicalErrorIf(not is_finite_vec(xs.sigma_t),
                       "CEPXS binary total cross section contains non-finite values.");
  OpenSnLogicalErrorIf(not is_finite_vec(xs.energy_deposition),
                       "CEPXS binary energy deposition contains non-finite values.");
  OpenSnLogicalErrorIf(not xs.stopping_power.empty() and not is_finite_vec(xs.stopping_power),
                       "CEPXS binary stopping power contains non-finite values.");
  OpenSnLogicalErrorIf(not xs.stopping_power.empty() and not IsNonNegative(xs.stopping_power),
                       "CEPXS binary stopping power contains negative values.");

  for (const auto& [g_begin, g_end] : particle_ranges)
  {
    OpenSnLogicalErrorIf(g_begin >= g_end or g_end >= xs.e_bounds.size(),
                         "CEPXS particle group range exceeds energy-bound storage.");
    for (size_t g = g_begin; g < g_end; ++g)
    {
      const double upper = g == g_begin ? xs.e_bounds.front() : xs.e_bounds[g];
      const double lower = xs.e_bounds[g + 1];
      const double pair_tol = 1.0e-12 * std::max(1.0, std::abs(upper));
      OpenSnLogicalErrorIf(upper <= lower + pair_tol,
                           "CEPXS particle groups must be strictly decreasing in energy. First "
                           "invalid pair in particle block at global index " +
                             std::to_string(g) + " -> " + std::to_string(g + 1) + " : " +
                             std::to_string(upper) + " <= " + std::to_string(lower) + ".");
    }
  }

  return xs;
}

} // namespace

MultiGroupXS
MultiGroupXS::LoadFromCEPXS(const std::string& filename,
                            int material_id,
                            bool csda_format,
                            const std::vector<std::string>& particle_order)
{
  MultiGroupXS mgxs;
  OpenSnLogicalErrorIf(not LooksLikeFortranBinary(filename),
                       "LoadFromCEPXS supports Fortran-record binary CEPXS only. File: \"" +
                         filename + "\".");
  const auto parsed = ParseCEPXSBFPBinary(
    filename, material_id, csda_format ? CEPXSRowFormat::CSDA : CEPXSRowFormat::LEGACY);

  mgxs.num_groups_ = parsed.num_groups;
  mgxs.scattering_order_ = parsed.scattering_order;
  mgxs.is_fissionable_ = parsed.is_fissionable;
  mgxs.num_precursors_ = mgxs.is_fissionable_ ? parsed.num_precursors : 0;

  mgxs.e_bounds_ = parsed.e_bounds;
  // BXSLIB stores one common upper energy followed by each group's lower energy.
  // At a species transition the preceding entry is the previous species' cutoff,
  // not the upper bound of the new species' first group.
  mgxs.e_upper_bounds_.reserve(mgxs.num_groups_);
  for (unsigned int g = 0; g < mgxs.num_groups_; ++g)
    mgxs.e_upper_bounds_.push_back(
      parsed.e_bounds[g] <= parsed.e_bounds[g + 1] ? parsed.e_bounds.front() : parsed.e_bounds[g]);
  mgxs.particle_types_ = ResolveParticleTypes(
    FindParticleGroupRanges(parsed.e_bounds), parsed.stopping_power, particle_order);
  mgxs.sigma_t_ = parsed.sigma_t;
  // Derive absorption from total and transfer matrices
  mgxs.sigma_a_.clear();
  mgxs.energy_deposition_ = parsed.energy_deposition;
  mgxs.custom_xs_["charge_deposition"] = parsed.charge_deposition;
  mgxs.custom_xs_["cepxs_charge_deposition"] = parsed.charge_deposition;
  mgxs.custom_xs_["cepxs_secondary_production"] = parsed.secondary_production;
  mgxs.stopping_power_ = parsed.stopping_power;
  mgxs.transfer_matrices_ = parsed.transfer_matrices;

  mgxs.ComputeAbsorption();
  mgxs.ComputeDiffusionParameters();
  return mgxs;
}

} // namespace opensn
