// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include "framework/utils/hdf_utils.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

using namespace opensn;

namespace
{

constexpr unsigned int kNumGroups = 2;

std::string
TempFileName(const std::string& stem, const std::string& extension)
{
  return (std::filesystem::temp_directory_path() /
          ("opensn_" + stem + "_" + std::to_string(getpid()) + "_r" +
           std::to_string(mpi_comm.rank()) + extension))
    .string();
}

std::shared_ptr<MultiGroupXS>
MakeTwoGroupXS()
{
  const auto path = TempFileName("uncollided_read_xs", ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS 2\nNUM_MOMENTS 2\n"
           "SIGMA_T_BEGIN\n0 1.0\n1 2.0\nSIGMA_T_END\n"
           "TRANSFER_MOMENTS_BEGIN\n"
           "M_GFROM_GTO_VAL 0 0 0 0.5\n"
           "M_GFROM_GTO_VAL 0 0 1 0.25\n"
           "M_GFROM_GTO_VAL 0 1 1 1.0\n"
           "M_GFROM_GTO_VAL 1 0 0 0.375\n"
           "M_GFROM_GTO_VAL 1 1 1 0.5\n"
           "TRANSFER_MOMENTS_END\n";
  }
  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path));
  std::filesystem::remove(path);
  return xs;
}

/// A 3D box of 4 x 4 x 4 hexahedra, two groups, P1 scattering (four moments), vacuum boundaries.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem()
{
  auto grid = BuildBoxMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 1}));
  groupset.AddParameter(
    "angular_quadrature",
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature3DXYZ>(2, 4, 1)));
  ParameterBlock groupsets("groupsets");
  groupsets.AddParameter(groupset);
  groupsets.ChangeToArray();

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeTwoGroupXS());
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", static_cast<int>(kNumGroups));
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  return DiscreteOrdinatesProblem::Create(params);
}

/// Exact (dyadic) uncollided flux of a cell, node, moment and group.
double
FluxValue(std::uint64_t global_id, size_t node, size_t moment, size_t group)
{
  return static_cast<double>(global_id) + 0.25 * static_cast<double>(node) +
         0.0078125 * static_cast<double>(group) + 1024.0 * static_cast<double>(moment);
}

/// The datasets of an uncollided flux file, one entry per cell in file order.
struct FileCells
{
  std::vector<std::uint64_t> ids;
  std::vector<std::uint64_t> node_counts;
  std::vector<std::vector<double>> x, y, z;
  std::vector<std::vector<double>> sigma_t;
};

/// The cells of all ranks, gathered on every rank in rank order. Collective.
FileCells
GatherCells(const DiscreteOrdinatesProblem& problem)
{
  const auto& grid = problem.GetGrid();
  const auto& discretization = problem.GetSpatialDiscretization();
  const auto& views = problem.GetCellTransportViews();
  std::vector<std::uint64_t> ids, counts;
  std::vector<double> x, y, z, sigma_t;
  for (const auto& cell : grid->GetLocalCells())
  {
    ids.push_back(cell->global_id);
    const auto num_nodes = discretization.GetCellNumNodes(*cell);
    counts.push_back(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto& vertex = grid->GlobalVertex(cell->vertex_ids[i]);
      x.push_back(vertex.x);
      y.push_back(vertex.y);
      z.push_back(vertex.z);
    }
    for (const auto value : views[cell->local_id].GetXS().GetSigmaTotal())
      sigma_t.push_back(value);
  }
  std::vector<std::uint64_t> all_ids, all_counts;
  std::vector<double> all_x, all_y, all_z, all_sigma_t;
  mpi_comm.all_gather(ids, all_ids);
  mpi_comm.all_gather(counts, all_counts);
  mpi_comm.all_gather(x, all_x);
  mpi_comm.all_gather(y, all_y);
  mpi_comm.all_gather(z, all_z);
  mpi_comm.all_gather(sigma_t, all_sigma_t);

  FileCells cells;
  size_t node = 0;
  for (size_t c = 0; c < all_ids.size(); ++c)
  {
    cells.ids.push_back(all_ids[c]);
    cells.node_counts.push_back(all_counts[c]);
    std::vector<double> cell_x, cell_y, cell_z, cell_sigma_t;
    for (size_t i = 0; i < all_counts[c]; ++i, ++node)
    {
      cell_x.push_back(all_x[node]);
      cell_y.push_back(all_y[node]);
      cell_z.push_back(all_z[node]);
    }
    for (size_t g = 0; g < kNumGroups; ++g)
      cell_sigma_t.push_back(all_sigma_t[c * kNumGroups + g]);
    cells.x.push_back(cell_x);
    cells.y.push_back(cell_y);
    cells.z.push_back(cell_z);
    cells.sigma_t.push_back(cell_sigma_t);
  }
  return cells;
}

/// Puts the cells in an order in which no two cells of a rank are likely to be neighbours in the
/// file: even global ids first, then odd ones, each in decreasing order.
void
Interleave(FileCells& cells)
{
  std::vector<size_t> order(cells.ids.size());
  for (size_t c = 0; c < order.size(); ++c)
    order[c] = c;
  std::sort(order.begin(),
            order.end(),
            [&cells](size_t a, size_t b)
            {
              const auto ka = cells.ids[a] % 2;
              const auto kb = cells.ids[b] % 2;
              return ka != kb ? ka < kb : cells.ids[a] > cells.ids[b];
            });
  FileCells sorted;
  for (const auto c : order)
  {
    sorted.ids.push_back(cells.ids[c]);
    sorted.node_counts.push_back(cells.node_counts[c]);
    sorted.x.push_back(cells.x[c]);
    sorted.y.push_back(cells.y[c]);
    sorted.z.push_back(cells.z[c]);
    sorted.sigma_t.push_back(cells.sigma_t[c]);
  }
  cells = sorted;
}

/// Writes an uncollided flux file of the cells, on rank 0, with the flux of FluxValue. Collective.
void
WriteFile(const std::string& file_name,
          const FileCells& cells,
          std::uint64_t global_cell_count,
          std::uint64_t moment_size_change = 0)
{
  if (mpi_comm.rank() == 0)
  {
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);
    std::vector<double> x, y, z, sigma_t;
    for (size_t c = 0; c < cells.ids.size(); ++c)
    {
      x.insert(x.end(), cells.x[c].begin(), cells.x[c].end());
      y.insert(y.end(), cells.y[c].begin(), cells.y[c].end());
      z.insert(z.end(), cells.z[c].begin(), cells.z[c].end());
      sigma_t.insert(sigma_t.end(), cells.sigma_t[c].begin(), cells.sigma_t[c].end());
    }
    ASSERT_TRUE(H5WriteDataset1D<std::uint64_t>(file.Id(), "cell ids", cells.ids));
    ASSERT_TRUE(H5WriteDataset1D<std::uint64_t>(file.Id(), "cell node counts", cells.node_counts));
    ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), "nodes x", x));
    ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), "nodes y", y));
    ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), "nodes z", z));
    ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), "cell sigma_t", sigma_t));
    const std::vector<std::string> moment_names = {"0,0", "1,-1", "1,0", "1,1"};
    for (size_t m = 0; m < moment_names.size(); ++m)
    {
      std::vector<double> values;
      for (size_t c = 0; c < cells.ids.size(); ++c)
        for (size_t i = 0; i < cells.node_counts[c]; ++i)
          for (size_t g = 0; g < kNumGroups; ++g)
            values.push_back(FluxValue(cells.ids[c], i, m, g));
      if (m == moment_names.size() - 1)
        values.resize(values.size() - moment_size_change);
      ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), moment_names[m], values));
    }
    ASSERT_TRUE(H5CreateAttribute<unsigned int>(file.Id(), "num groups", kNumGroups));
    ASSERT_TRUE(H5CreateAttribute<unsigned int>(file.Id(), "max moment order", 1U));
    ASSERT_TRUE(
      H5CreateAttribute<std::uint64_t>(file.Id(), "global cell count", global_cell_count));
    ASSERT_TRUE(H5CreateAttribute<double>(file.Id(), "production", 1.0));
    ASSERT_TRUE(H5CreateAttribute<double>(file.Id(), "out-flow", 0.5));
    ASSERT_TRUE(H5CreateAttribute<std::uint64_t>(file.Id(), "reflecting boundary count", 0));
  }
  mpi_comm.barrier();
}

/// The file name every rank uses: the name rank 0 writes.
std::string
SharedFileName(const std::string& stem)
{
  int pid = getpid();
  mpi_comm.broadcast(pid, 0);
  return (std::filesystem::temp_directory_path() /
          ("opensn_" + stem + "_" + std::to_string(pid) + ".h5"))
    .string();
}

} // namespace

/// Every rank reads its own cells, at the positions of the problem's unknowns, from a file whose
/// cells are in an order unrelated to the partition, so a rank's cells are scattered over the
/// file.
TEST(UncollidedFluxReadTest, ReadsLocalCellsInAnyFileOrder)
{
  const auto problem = MakeProblem();
  auto cells = GatherCells(*problem);
  Interleave(cells);
  const auto global_cell_count = problem->GetGrid()->GetGlobalNumberOfCells();
  const auto file_name = SharedFileName("uncollided_read");
  WriteFile(file_name, cells, global_cell_count);

  const auto data = DiscreteOrdinatesProblemIO::ReadUncollidedFlux(*problem, file_name);
  EXPECT_EQ(data.num_groups, kNumGroups);
  EXPECT_EQ(data.max_moment_order, 1U);
  EXPECT_EQ(data.global_cell_count, global_cell_count);
  EXPECT_EQ(data.source_rate, 1.0);
  EXPECT_EQ(data.outflow_rate, 0.5);

  const auto& discretization = problem->GetSpatialDiscretization();
  const auto& views = problem->GetCellTransportViews();
  ASSERT_EQ(data.local_flux_moments.size(), problem->GetPhiNewLocal().size());
  size_t mismatches = 0;
  for (const auto& cell : problem->GetGrid()->GetLocalCells())
    for (size_t i = 0; i < discretization.GetCellNumNodes(*cell); ++i)
      for (size_t m = 0; m < 4; ++m)
        for (size_t g = 0; g < kNumGroups; ++g)
          if (data.local_flux_moments[views[cell->local_id].MapDOF(i, m, g)] !=
              FluxValue(cell->global_id, i, m, g))
            ++mismatches;
  EXPECT_EQ(mismatches, 0U);

  mpi_comm.barrier();
  if (mpi_comm.rank() == 0)
    std::filesystem::remove(file_name);
}

/// Files that do not describe the mesh are rejected: a missing cell, a cell listed twice, a moved
/// node, and a moment dataset of the wrong size. Each check runs on the rank that owns the cell,
/// so these cases run on one rank.
TEST(UncollidedFluxReadTest, RejectsInconsistentFiles)
{
  if (mpi_comm.size() != 1)
    GTEST_SKIP() << "The rejections are rank-local, so they are tested on one rank.";

  const auto problem = MakeProblem();
  const auto cells = GatherCells(*problem);
  const auto global_cell_count = problem->GetGrid()->GetGlobalNumberOfCells();
  const auto file_name = TempFileName("uncollided_reject", ".h5");

  const auto expect_rejected =
    [&](const FileCells& file_cells, std::uint64_t moment_size_change, const std::string& expected)
  {
    WriteFile(file_name, file_cells, global_cell_count, moment_size_change);
    try
    {
      DiscreteOrdinatesProblemIO::ReadUncollidedFlux(*problem, file_name);
      ADD_FAILURE() << "accepted a file that should fail with \"" << expected << "\"";
    }
    catch (const std::invalid_argument& error)
    {
      EXPECT_NE(std::string(error.what()).find(expected), std::string::npos) << error.what();
    }
  };

  // One cell fewer than the mesh.
  auto missing = cells;
  for (auto* entries : {&missing.x, &missing.y, &missing.z, &missing.sigma_t})
    entries->pop_back();
  missing.ids.pop_back();
  missing.node_counts.pop_back();
  expect_rejected(missing, 0, "cells, but its mesh has");

  // The last cell replaced by a second copy of the first.
  auto duplicate = cells;
  duplicate.ids.back() = duplicate.ids.front();
  expect_rejected(duplicate, 0, "duplicate cell id");

  // A node of one cell moved.
  auto moved = cells;
  moved.x[3][1] += 1.0e-6;
  expect_rejected(moved, 0, "mesh coordinate mismatch");

  // The last moment dataset one value short.
  expect_rejected(cells, 1, "has size");

  std::filesystem::remove(file_name);
}
