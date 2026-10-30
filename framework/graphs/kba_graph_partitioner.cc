// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <sstream>

namespace opensn
{

InputParameters
KBAGraphPartitioner::GetInputParameters()
{
  InputParameters params = GraphPartitioner::GetInputParameters();

  params.AddOptionalParameter("nx", 1, "The number of partitions in x");
  params.AddOptionalParameter("ny", 1, "The number of partitions in y");
  params.AddOptionalParameter("nz", 1, "The number of partitions in z");

  params.AddOptionalParameter(
    "xcuts", std::vector<double>{}, "Location of the internal x-cuts. Require nx-1 entries");
  params.AddOptionalParameter(
    "ycuts", std::vector<double>{}, "Location of the internal y-cuts. Require ny-1 entries");
  params.AddOptionalParameter(
    "zcuts", std::vector<double>{}, "Location of the internal z-cuts. Require nz-1 entries");

  params.ConstrainParameterRange("nx", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("ny", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("nz", AllowableRangeLowLimit::New(1));

  return params;
}

std::shared_ptr<KBAGraphPartitioner>
KBAGraphPartitioner::Create(const ParameterBlock& params)
{
  return CreateObject<KBAGraphPartitioner>("mesh::KBAGraphPartitioner", params);
}

KBAGraphPartitioner::KBAGraphPartitioner(const InputParameters& params)
  : GraphPartitioner(params),
    nx_(params.GetParamValue<size_t>("nx")),
    ny_(params.GetParamValue<size_t>("ny")),
    nz_(params.GetParamValue<size_t>("nz")),
    xcuts_(params.GetParamVectorValue<double>("xcuts")),
    ycuts_(params.GetParamVectorValue<double>("ycuts")),
    zcuts_(params.GetParamVectorValue<double>("zcuts")),
    coordinate_infos_{CoordinateInfo{&xcuts_, nx_, "x"},
                      CoordinateInfo{&ycuts_, ny_, "y"},
                      CoordinateInfo{&zcuts_, nz_, "z"}}
{
  for (const auto& [cuts_ptr, n, name] : coordinate_infos_)
  {
    const auto& cuts = *cuts_ptr;

    // Check number of items
    if (cuts.size() != (n - 1))
    {
      std::ostringstream err;
      err << "The number of cuts supplied for \"" << name << "cuts\" is not equal to n" << name
          << "-1.";
      OpenSnInvalidArgument(err.str());
    }
    if (cuts.empty())
      continue;

    // Cuts must be finite and strictly increasing
    for (size_t c = 0; c < cuts.size(); ++c)
      OpenSnInvalidArgumentIf(not std::isfinite(cuts[c]) or (c > 0 and cuts[c] <= cuts[c - 1]),
                              "KBAGraphPartitioner: " + std::string(name) +
                                "cuts must be finite and strictly increasing.");
  } // for each coordinate
}

namespace
{

/**
 * Size, in cellsets, of equal sub-grids of the n[0] x n[1] x n[2] cellset array that give each of
 * num_parts partitions a brick-shaped subdomain of whole cellsets, or nullopt when there is none.
 * Consecutive cellset ids are preferred when they already form such sub-grids. Otherwise, the
 * sub-grid with the least surface is used.
 */
std::optional<std::array<size_t, 3>>
FindCellsetGrouping(const std::array<size_t, 3>& n, size_t num_parts)
{
  const size_t num_cellsets = n[0] * n[1] * n[2];
  if (num_parts == 0 or num_cellsets < num_parts or num_cellsets % num_parts != 0)
    return std::nullopt;
  const size_t per_part = num_cellsets / num_parts;

  // Consecutive cellset ids (x fastest) form a brick when they fill part of a row, whole rows of a
  // layer, or whole layers.
  if (n[0] % per_part == 0)
    return std::array<size_t, 3>{per_part, 1, 1};
  if (per_part % n[0] == 0 and n[1] % (per_part / n[0]) == 0)
    return std::array<size_t, 3>{n[0], per_part / n[0], 1};
  if (per_part % (n[0] * n[1]) == 0 and n[2] % (per_part / (n[0] * n[1])) == 0)
    return std::array<size_t, 3>{n[0], n[1], per_part / (n[0] * n[1])};

  std::optional<std::array<size_t, 3>> best;
  size_t best_surface = std::numeric_limits<size_t>::max();
  for (size_t gx = n[0]; gx > 0; --gx)
  {
    if (n[0] % gx != 0 or per_part % gx != 0)
      continue;
    for (size_t gy = n[1]; gy > 0; --gy)
    {
      if (n[1] % gy != 0 or (per_part / gx) % gy != 0)
        continue;
      const size_t gz = per_part / (gx * gy);
      if (n[2] % gz != 0)
        continue;
      const size_t surface = gx * gy + gy * gz + gz * gx;
      if (surface < best_surface)
      {
        best_surface = surface;
        best = std::array<size_t, 3>{gx, gy, gz};
      }
    }
  }
  return best;
}

} // namespace

std::vector<int>
KBAGraphPartitioner::Partition(const std::vector<std::vector<uint64_t>>& graph,
                               const std::vector<Vector3>& centroids,
                               int number_of_parts)
{
  log.Log0Verbose1() << "Partitioning with KBAGraphPartitioner";

  OpenSnLogicalErrorIf(centroids.size() != graph.size(),
                       "Graph number of entries not equal to centroids' number of entries.");
  OpenSnInvalidArgumentIf(number_of_parts < 1,
                          "KBAGraphPartitioner: the number of partitions must be at least 1.");
  const auto num_parts = static_cast<size_t>(number_of_parts);
  const size_t num_cells = graph.size();

  // Cutting along an axis on which every centroid has the same coordinate (z for a 2D mesh, x and
  // y for a 1D mesh) leaves every slab but one without cells. Centroids of cells in one plane can
  // differ in the last bits when the plane is not at zero, so the comparison allows a few units in
  // the last place of the coordinate itself. Any real extent, however thin, is larger.
  if (num_cells > 0)
  {
    for (size_t d = 0; d < 3; ++d)
    {
      const auto [min_it, max_it] =
        std::minmax_element(centroids.begin(),
                            centroids.end(),
                            [d](const Vector3& a, const Vector3& b) { return a[d] < b[d]; });
      const double lo = (*min_it)[d];
      const double hi = (*max_it)[d];
      const double rounding =
        64.0 * std::numeric_limits<double>::epsilon() * std::max(std::fabs(lo), std::fabs(hi));
      const auto& [cuts, n, name] = coordinate_infos_[d];
      if (n > 1 and hi - lo <= rounding)
      {
        std::ostringstream err;
        err << "KBAGraphPartitioner: n" << name << " = " << n
            << ", but every cell centroid has the same " << name << " coordinate (for example, "
            << (d == 2 ? "a 2D mesh" : "a 1D mesh") << "), so the cells cannot be partitioned in "
            << name << ". Use n" << name << " = 1 and no " << name << "cuts.";
        OpenSnInvalidArgument(err.str());
      }
    }
  }

  std::vector<size_t> pids(num_cells, 0);
  std::array<std::vector<size_t>, 3> cells_per_slab;
  for (size_t d = 0; d < 3; ++d)
    cells_per_slab[d].assign(coordinate_infos_[d].n_, 0);
  for (size_t c = 0; c < num_cells; ++c)
  {
    const auto& point = centroids[c];
    // Partitions the point
    std::array<size_t, 3> p_vals = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i)
    {
      const auto& cuts = *coordinate_infos_[i].cuts_;
      const size_t num_cuts = cuts.size();

      size_t p_val = 0;
      bool home_found = false;
      for (size_t j = 0; j < num_cuts; ++j)
        if (cuts[j] > point[i])
        {
          p_val = j;
          home_found = true;
          break;
        }

      p_vals[i] = home_found ? p_val : (coordinate_infos_[i].n_ - 1);
    }

    const auto nx = static_cast<int64_t>(coordinate_infos_[0].n_);
    const auto ny = static_cast<int64_t>(coordinate_infos_[1].n_);

    const auto i = static_cast<int64_t>(p_vals[0]);
    const auto j = static_cast<int64_t>(p_vals[1]);
    const auto k = static_cast<int64_t>(p_vals[2]);

    pids[c] = nx * ny * k + nx * j + i;
    for (size_t d = 0; d < 3; ++d)
      ++cells_per_slab[d][p_vals[d]];
  } // for cell c

  // A slab without cells means a cut lies outside the mesh or two cuts have no cell centroid
  // between them. Its cellsets would be empty, and refilling them breaks the brick decomposition.
  for (size_t d = 0; d < 3; ++d)
  {
    const auto& [cuts, n, name] = coordinate_infos_[d];
    for (size_t slab = 0; slab < n; ++slab)
      if (num_cells > 0 and cells_per_slab[d][slab] == 0)
      {
        std::ostringstream err;
        err << "KBAGraphPartitioner: no cell centroid lies in " << name << " slab " << slab
            << " of " << n << ". Every slab must contain cells. Check " << name
            << "cuts against the mesh extent.";
        OpenSnInvalidArgument(err.str());
      }
  }

  const std::array<size_t, 3> num_cellsets_per_axis = {nx_, ny_, nz_};
  const size_t num_cellsets = nx_ * ny_ * nz_;
  std::vector<int> real_pids(num_cells, 0);
  if (const auto group = FindCellsetGrouping(num_cellsets_per_axis, num_parts))
  {
    // Each partition owns one gx x gy x gz sub-grid of cellsets.
    const auto [gx, gy, gz] = *group;
    const size_t groups_x = nx_ / gx;
    const size_t groups_y = ny_ / gy;
    for (size_t c = 0; c < num_cells; ++c)
    {
      const size_t i = pids[c] % nx_;
      const size_t j = (pids[c] / nx_) % ny_;
      const size_t k = pids[c] / (nx_ * ny_);
      real_pids[c] = static_cast<int>(i / gx + groups_x * (j / gy + groups_y * (k / gz)));
    }

    // Cellsets can be empty when the mesh does not fill the overlay (e.g. an L-shaped domain).
    std::vector<size_t> cells_per_part(num_parts, 0);
    for (const int pid : real_pids)
      ++cells_per_part[pid];
    const auto num_empty = std::count(cells_per_part.begin(), cells_per_part.end(), size_t{0});
    if (num_empty > 0 and num_cells >= num_parts)
      log.Log0Warning() << "KBAGraphPartitioner: " << num_empty
                        << " partition(s) own only cellsets without cells. They take cells from "
                           "the largest partitions, whose subdomains are then no longer bricks.";
  }
  else
  {
    if (num_cellsets < num_parts)
      log.Log0Warning() << "KBAGraphPartitioner: nx*ny*nz = " << num_cellsets << " cellsets for "
                        << number_of_parts
                        << " partitions. The empty partitions take cells from the largest ones.";
    else
      log.Log0Warning() << "KBAGraphPartitioner: the nx*ny*nz = " << num_cellsets
                        << " cellsets cannot be grouped into " << number_of_parts
                        << " equal brick-shaped subdomains, so partitions take consecutive "
                           "cellsets, which may not form a brick and can create cyclic sweep "
                           "dependencies.";

    const auto pid_subsets = MakeSubSets(num_cellsets, num_parts);
    for (size_t c = 0; c < num_cells; ++c)
    {
      for (int p = 0; p < number_of_parts; ++p)
      {
        if (pids[c] >= pid_subsets[p].ss_begin and pids[c] <= pid_subsets[p].ss_end)
          real_pids[c] = p;
      }
    }
  }

  log.Log0Verbose1() << "Done partitioning with KBAGraphPartitioner";

  return real_pids;
}

} // namespace opensn
