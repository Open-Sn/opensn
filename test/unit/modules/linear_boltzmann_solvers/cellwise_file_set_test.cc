// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/cellwise_file_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/mpi/record_exchange.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include "framework/utils/hdf_utils.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

std::shared_ptr<MultiGroupXS>
MakeTwoGroupXS(const std::filesystem::path& directory)
{
  const auto path = directory / ("xs_r" + std::to_string(mpi_comm.rank()) + ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS 2\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n0 1.0\n1 2.0\nSIGMA_T_END\n";
  }
  return std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
}

/// A 2D XY problem on an 8 x 8 mesh with two groups, one groupset per group, and P1 moments.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem(const std::filesystem::path& directory)
{
  auto grid = BuildSquareMesh(1.0, 8, 0.0);
  grid->SetUniformBlockID(0);

  const auto quadrature =
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature2DXY>(2, 4, 1));
  ParameterBlock groupsets("groupsets");
  for (int g = 0; g < 2; ++g)
  {
    ParameterBlock groupset(std::to_string(g));
    groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{g, g}));
    groupset.AddParameter("angular_quadrature", quadrature);
    groupsets.AddParameter(groupset);
  }
  groupsets.ChangeToArray();

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeTwoGroupXS(directory));
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 2);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  return DiscreteOrdinatesProblem::Create(params);
}

/// A cell of the mesh and its node locations.
struct CellNodes
{
  std::uint64_t id = 0;
  std::vector<Vector3> nodes;
};

/// The cells of the whole mesh, gathered on rank 0 in global id order; empty on other ranks.
std::vector<CellNodes>
GatherCells(const DiscreteOrdinatesProblem& problem)
{
  const auto& discretization = problem.GetSpatialDiscretization();
  RecordBatch batch(1);
  for (const auto& cell : problem.GetGrid()->GetLocalCells())
  {
    batch.AddRecord(0)[0] = static_cast<std::int64_t>(cell->global_id);
    for (const auto& node : discretization.GetCellNodeLocations(*cell))
    {
      const double xyz[3] = {node.x, node.y, node.z};
      batch.AppendPayload(xyz, 3);
    }
  }
  const auto received = ExchangeRecords(mpi_comm, std::move(batch), "Gathering the mesh");
  std::vector<CellNodes> cells;
  for (std::size_t r = 0; r < received.GetNumRecords(); ++r)
  {
    auto& cell = cells.emplace_back();
    cell.id = static_cast<std::uint64_t>(received.GetHeader(r)[0]);
    const double* xyz = received.GetPayload(r);
    for (std::size_t i = 0; i < received.GetPayloadSize(r) / 3; ++i)
      cell.nodes.emplace_back(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);
  }
  std::sort(cells.begin(), cells.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
  return cells;
}

/**
 * Value of unknown `a` (a moment or a direction) and group `g` of groupset `gs` at `node`. Node
 * coordinates are multiples of 1/8, so every value is exact and survives copying unchanged.
 */
double
Value(const Vector3& node, const std::size_t a, const std::size_t g, const std::size_t gs)
{
  return node.x + 2.0 * node.y + 4.0 * static_cast<double>(a) + 256.0 * static_cast<double>(g) +
         1024.0 * static_cast<double>(gs);
}

/// How a synthetic file set is written.
struct SetSpec
{
  bool angular = false;
  std::uint64_t num_files = 1;
  bool single_file = false;
  /// The num_partitions attribute of each file, or nullopt to omit it; the default is num_files.
  std::function<std::optional<std::uint64_t>(std::uint64_t)> recorded_files;
  /// The files that hold each cell; the default is file (7 id) mod num_files.
  std::function<std::vector<std::uint64_t>(std::uint64_t)> files_of;
  /// Applied to the nodes of each written cell.
  std::function<void(std::uint64_t, std::vector<Vector3>&)> edit_nodes;
  /// Applied to each written value.
  std::function<double(std::uint64_t, double)> edit_value;
  /// A cell id that is not in the mesh, written to file 0 with the nodes of cell 0, or 0 for none.
  std::uint64_t unknown_cell = 0;
  unsigned int role = static_cast<unsigned int>(FluxSolutionRole::FORWARD);
  unsigned int orientation = static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED);
  /// The (role, orientation) of each file, overriding the two above, or nullopt to omit both.
  std::function<std::optional<std::pair<unsigned int, unsigned int>>(std::uint64_t)> provenance_of;
  /// Omit the quadrature weights of each groupset of an angular set.
  bool omit_weights = false;
};

class CellwiseFileSetTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    std::uint64_t nonce =
      static_cast<std::uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count());
    mpi_comm.broadcast(nonce, 0);
    directory_ = std::filesystem::temp_directory_path() /
                 ("opensn_cellwise_file_set_" + std::to_string(nonce));
    if (mpi_comm.rank() == 0)
      std::filesystem::create_directories(directory_);
    mpi_comm.barrier();
    problem_ = MakeProblem(directory_);
  }

  void TearDown() override
  {
    mpi_comm.barrier();
    if (mpi_comm.rank() == 0)
      std::filesystem::remove_all(directory_);
  }

  /// Write a file set described by `spec` on rank 0 and return its base. Collective.
  std::string Write(const std::string& name, const SetSpec& spec)
  {
    const auto base = (directory_ / name).string();
    const auto cells = GatherCells(*problem_);
    if (mpi_comm.rank() == 0)
      for (std::uint64_t file_id = 0; file_id < spec.num_files; ++file_id)
        WriteFile(spec, base, file_id, cells);
    mpi_comm.barrier();
    return base;
  }

  /// Expect `phi` to hold the flux moments of every local node, exactly.
  void ExpectFluxMoments(const std::vector<double>& phi) const
  {
    const auto& discretization = problem_->GetSpatialDiscretization();
    const auto& uk_man = problem_->GetUnknownManager();
    ASSERT_EQ(phi.size(), discretization.GetNumLocalDOFs(uk_man));
    for (const auto& cell : problem_->GetGrid()->GetLocalCells())
    {
      const auto& nodes = discretization.GetCellNodeLocations(*cell);
      for (std::size_t i = 0; i < nodes.size(); ++i)
        for (unsigned int m = 0; m < problem_->GetNumMoments(); ++m)
          for (unsigned int g = 0; g < 2; ++g)
            ASSERT_EQ(phi[discretization.MapDOFLocal(*cell, i, uk_man, m, g)],
                      Value(nodes[i], m, g, 0))
              << "cell " << cell->global_id << " node " << i;
    }
  }

  /// Expect `psi` to hold the angular fluxes of every groupset and local node, exactly.
  void ExpectAngularFluxes(const std::vector<std::vector<double>>& psi) const
  {
    const auto& discretization = problem_->GetSpatialDiscretization();
    const auto& groupsets = problem_->GetGroupsets();
    ASSERT_EQ(psi.size(), groupsets.size());
    for (std::size_t gs = 0; gs < groupsets.size(); ++gs)
    {
      const auto& uk_man = groupsets[gs].psi_uk_man_;
      ASSERT_EQ(psi[gs].size(), discretization.GetNumLocalDOFs(uk_man));
      for (const auto& cell : problem_->GetGrid()->GetLocalCells())
      {
        const auto& nodes = discretization.GetCellNodeLocations(*cell);
        for (std::size_t i = 0; i < nodes.size(); ++i)
          for (std::size_t n = 0; n < groupsets[gs].quadrature->GetNumAngles(); ++n)
            ASSERT_EQ(psi[gs][discretization.MapDOFLocal(*cell, i, uk_man, n, 0)],
                      Value(nodes[i], n, 0, gs))
              << "groupset " << gs << " cell " << cell->global_id << " node " << i;
      }
    }
  }

  std::filesystem::path directory_;
  std::shared_ptr<DiscreteOrdinatesProblem> problem_;

private:
  void WriteFile(const SetSpec& spec,
                 const std::string& base,
                 const std::uint64_t file_id,
                 const std::vector<CellNodes>& cells) const
  {
    // The cells of this file, nodes reversed in every third cell, so node order differs.
    std::vector<CellNodes> file_cells;
    for (const auto& cell : cells)
    {
      const auto files = spec.files_of ? spec.files_of(cell.id)
                                       : std::vector<std::uint64_t>{(7 * cell.id) % spec.num_files};
      if (std::find(files.begin(), files.end(), file_id) == files.end())
        continue;
      auto& file_cell = file_cells.emplace_back(cell);
      if (cell.id % 3 == 0)
        std::reverse(file_cell.nodes.begin(), file_cell.nodes.end());
      if (spec.edit_nodes)
        spec.edit_nodes(cell.id, file_cell.nodes);
    }
    if (spec.unknown_cell != 0 and file_id == 0)
      file_cells.push_back({spec.unknown_cell, cells.front().nodes});

    std::vector<std::uint64_t> cell_ids, num_cell_nodes;
    std::vector<double> nodes_x, nodes_y, nodes_z;
    for (const auto& cell : file_cells)
    {
      cell_ids.push_back(cell.id);
      num_cell_nodes.push_back(cell.nodes.size());
      for (const auto& node : cell.nodes)
      {
        nodes_x.push_back(node.x);
        nodes_y.push_back(node.y);
        nodes_z.push_back(node.z);
      }
    }
    // Values in file order: cell, node, unknown, group.
    const auto Values =
      [&](const std::size_t num_unknowns, const std::size_t num_groups, const std::size_t gs)
    {
      std::vector<double> values;
      for (const auto& cell : file_cells)
        for (const auto& node : cell.nodes)
          for (std::size_t a = 0; a < num_unknowns; ++a)
            for (std::size_t g = 0; g < num_groups; ++g)
            {
              const auto value = Value(node, a, g, gs);
              values.push_back(spec.edit_value ? spec.edit_value(cell.id, value) : value);
            }
      return values;
    };

    const auto file_name = spec.single_file ? base + ".h5" : PartitionFileName(base, file_id);
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);
    const auto id = file.Id();
    const auto recorded =
      spec.recorded_files ? spec.recorded_files(file_id) : std::optional(spec.num_files);
    if (recorded.has_value())
    {
      ASSERT_TRUE(H5CreateAttribute(id, "num_partitions", *recorded));
    }
    const auto provenance = spec.provenance_of
                              ? spec.provenance_of(file_id)
                              : std::optional(std::make_pair(spec.role, spec.orientation));
    if (provenance.has_value())
    {
      ASSERT_TRUE(H5CreateAttribute(id, "solution_role", provenance->first));
      ASSERT_TRUE(H5CreateAttribute(id, "adjoint_orientation", provenance->second));
    }
    ASSERT_TRUE(H5CreateGroup(id, "mesh"));
    ASSERT_TRUE(H5CreateAttribute(id, "mesh/num_local_cells", std::uint64_t{cell_ids.size()}));
    ASSERT_TRUE(H5CreateAttribute(id, "mesh/num_local_nodes", std::uint64_t{nodes_x.size()}));
    ASSERT_TRUE(H5WriteDataset1D(id, "mesh/cell_ids", cell_ids));
    ASSERT_TRUE(H5WriteDataset1D(id, "mesh/num_cell_nodes", num_cell_nodes));
    ASSERT_TRUE(H5WriteDataset1D(id, "mesh/nodes_x", nodes_x));
    ASSERT_TRUE(H5WriteDataset1D(id, "mesh/nodes_y", nodes_y));
    ASSERT_TRUE(H5WriteDataset1D(id, "mesh/nodes_z", nodes_z));

    if (not spec.angular)
    {
      ASSERT_TRUE(H5CreateAttribute(id, "num_moments", problem_->GetNumMoments()));
      ASSERT_TRUE(H5CreateAttribute(id, "num_groups", 2U));
      ASSERT_TRUE(H5WriteDataset1D(id, "values", Values(problem_->GetNumMoments(), 2, 0)));
      return;
    }
    const auto& groupsets = problem_->GetGroupsets();
    ASSERT_TRUE(H5CreateAttribute(id, "num_groupsets", std::uint64_t{groupsets.size()}));
    for (std::size_t gs = 0; gs < groupsets.size(); ++gs)
    {
      const auto& quadrature = groupsets[gs].quadrature;
      const auto group = "groupset_" + std::to_string(groupsets[gs].id);
      std::vector<double> omegas;
      for (const auto& omega : quadrature->GetOmegas())
        omegas.insert(omegas.end(), {omega.x, omega.y, omega.z});
      ASSERT_TRUE(H5CreateGroup(id, group));
      ASSERT_TRUE(H5CreateAttribute(
        id, group + "/num_directions", std::uint64_t{quadrature->GetNumAngles()}));
      ASSERT_TRUE(H5CreateAttribute(id, group + "/num_groups", 1U));
      ASSERT_TRUE(H5WriteDataset1D(id, group + "/omegas", omegas));
      if (not spec.omit_weights)
      {
        ASSERT_TRUE(H5WriteDataset1D(id, group + "/weights", quadrature->GetWeights()));
      }
      ASSERT_TRUE(
        H5WriteDataset1D(id, group + "/values", Values(quadrature->GetNumAngles(), 1, gs)));
    }
  }
};

} // namespace

TEST_F(CellwiseFileSetTest, ReadsFluxMomentsWrittenByAnyRankCount)
{
  for (const std::uint64_t num_files : {1, 2, 3, 5})
  {
    SetSpec spec;
    spec.num_files = num_files;
    const auto base = Write("phi" + std::to_string(num_files) + "_", spec);
    std::vector<double> phi;
    LBSSolverIO::ReadFluxMoments(*problem_, base, false, phi);
    ExpectFluxMoments(phi);
  }

  // The problem's own flux moments are the default destination.
  SetSpec spec;
  spec.num_files = 3;
  LBSSolverIO::ReadFluxMoments(*problem_, Write("phi_default_", spec), false);
  ExpectFluxMoments(problem_->GetPhiOldLocal());
}

TEST_F(CellwiseFileSetTest, ReadsAngularFluxesWrittenByAnyRankCount)
{
  for (const std::uint64_t num_files : {1, 2, 5})
  {
    SetSpec spec;
    spec.angular = true;
    spec.num_files = num_files;
    const auto base = Write("psi" + std::to_string(num_files) + "_", spec);
    std::vector<std::vector<double>> psi;
    DiscreteOrdinatesProblemIO::ReadAngularFluxes(*problem_, base, psi);
    ExpectAngularFluxes(psi);
  }

  SetSpec spec;
  spec.angular = true;
  spec.num_files = 2;
  DiscreteOrdinatesProblemIO::ReadAngularFluxes(*problem_, Write("psi_default_", spec));
  ExpectAngularFluxes(problem_->GetPsiNewLocal());
}

TEST_F(CellwiseFileSetTest, SmallSendBudgetGivesTheSameValues)
{
  SetSpec spec;
  spec.num_files = 3;
  const auto base = Write("phi_", spec);
  const auto& discretization = problem_->GetSpatialDiscretization();
  const auto& uk_man = problem_->GetUnknownManager();
  const auto num_moments = problem_->GetNumMoments();
  // A one-byte budget sends one cell per rank per exchange.
  std::vector<double> phi(discretization.GetNumLocalDOFs(uk_man), 0.0);
  ReadCellwiseFileSet(
    base,
    false,
    {{"values", 2ULL * num_moments}},
    *problem_->GetGrid(),
    discretization,
    [](hid_t, const std::string&) {},
    [&](std::size_t, const Cell& cell, const std::size_t node, const double* values)
    {
      for (unsigned int m = 0; m < num_moments; ++m)
        for (unsigned int g = 0; g < 2; ++g)
          phi[discretization.MapDOFLocal(cell, node, uk_man, m, g)] = values[2 * m + g];
    },
    mpi_comm,
    "Reading test flux moments",
    1);
  ExpectFluxMoments(phi);
}

TEST_F(CellwiseFileSetTest, ReadsWrittenFilesOnTheirOwnPartition)
{
  // A write followed by a read on the same ranks gives back the same vector.
  auto& phi = problem_->GetPhiNewLocal();
  for (std::size_t k = 0; k < phi.size(); ++k)
    phi[k] = 0.25 * static_cast<double>(k % 97);
  const auto base = (directory_ / "written_").string();
  LBSSolverIO::WriteFluxMoments(*problem_, base);
  mpi_comm.barrier();
  std::vector<double> read;
  LBSSolverIO::ReadFluxMoments(*problem_, base, false, read);
  EXPECT_EQ(read, phi);
}

TEST_F(CellwiseFileSetTest, ReadsSetsWithoutAFileCountOnTheirRankCount)
{
  // Files of earlier versions do not record the number of files; a set written by as many ranks
  // as read it is still read.
  SetSpec spec;
  spec.num_files = static_cast<std::uint64_t>(mpi_comm.size());
  spec.recorded_files = [](std::uint64_t) { return std::nullopt; };
  std::vector<double> phi;
  LBSSolverIO::ReadFluxMoments(*problem_, Write("unrecorded_", spec), false, phi);
  ExpectFluxMoments(phi);

  // Written by more ranks, its last files are not read, and their cells are missing. Cell 0 is in
  // the last file, which is not read, whatever the rank count.
  spec.num_files = static_cast<std::uint64_t>(mpi_comm.size()) + 1;
  spec.files_of = [last = spec.num_files - 1](std::uint64_t id)
  { return std::vector<std::uint64_t>{id == 0 ? last : 0}; };
  EXPECT_THROW(LBSSolverIO::ReadFluxMoments(*problem_, Write("unrecorded_more_", spec), false, phi),
               std::runtime_error);
}

TEST_F(CellwiseFileSetTest, ReadsASingleFileOnAnyRankCount)
{
  SetSpec spec;
  spec.single_file = true;
  std::vector<double> phi;
  LBSSolverIO::ReadFluxMoments(*problem_, Write("single", spec), true, phi);
  ExpectFluxMoments(phi);

  spec.recorded_files = [](std::uint64_t) { return std::nullopt; };
  LBSSolverIO::ReadFluxMoments(*problem_, Write("single_unrecorded", spec), true, phi);
  ExpectFluxMoments(phi);

  // A file of a larger set is not a single file.
  spec.recorded_files = [](std::uint64_t) { return std::optional<std::uint64_t>(3); };
  EXPECT_THROW(LBSSolverIO::ReadFluxMoments(*problem_, Write("single_of_3", spec), true, phi),
               std::runtime_error);
}

TEST_F(CellwiseFileSetTest, RejectsInconsistentSetsAndLeavesTheDestination)
{
  const std::vector<double> untouched(4, -1.0);
  const auto ExpectRejected = [&](const std::string& name, const SetSpec& spec)
  {
    const auto base = spec.num_files == 0 ? (directory_ / name).string() : Write(name, spec);
    auto phi = untouched;
    EXPECT_THROW(LBSSolverIO::ReadFluxMoments(*problem_, base, false, phi), std::runtime_error)
      << name;
    EXPECT_EQ(phi, untouched) << name;
  };

  SetSpec missing;
  missing.num_files = 0;
  ExpectRejected("missing_", missing);

  SetSpec counts;
  counts.num_files = 2;
  counts.recorded_files = [](std::uint64_t file) { return std::optional<std::uint64_t>(file + 2); };
  ExpectRejected("counts_", counts);

  SetSpec mixed;
  mixed.num_files = 2;
  mixed.recorded_files = [](std::uint64_t file)
  { return file == 0 ? std::optional<std::uint64_t>(2) : std::nullopt; };
  ExpectRejected("mixed_", mixed);

  SetSpec duplicate;
  duplicate.num_files = 2;
  duplicate.files_of = [](std::uint64_t id)
  { return id == 5 ? std::vector<std::uint64_t>{0, 1} : std::vector<std::uint64_t>{id % 2}; };
  ExpectRejected("duplicate_", duplicate);

  SetSpec absent;
  absent.num_files = 2;
  absent.files_of = [](std::uint64_t id)
  { return id == 5 ? std::vector<std::uint64_t>{} : std::vector<std::uint64_t>{id % 2}; };
  ExpectRejected("absent_", absent);

  SetSpec unknown;
  unknown.num_files = 2;
  unknown.unknown_cell = 1000000;
  ExpectRejected("unknown_", unknown);

  SetSpec moved;
  moved.num_files = 2;
  moved.edit_nodes = [](std::uint64_t id, std::vector<Vector3>& nodes)
  {
    if (id == 5)
      nodes[0].x += 1.0e-3;
  };
  ExpectRejected("moved_", moved);

  SetSpec non_finite;
  non_finite.num_files = 2;
  non_finite.edit_value = [](std::uint64_t id, double value)
  { return id == 5 ? std::numeric_limits<double>::quiet_NaN() : value; };
  ExpectRejected("non_finite_", non_finite);

  // An infinite coordinate would make the node-matching tolerance infinite. The values, which are
  // written from the node locations, are kept finite so that only the coordinate is invalid.
  SetSpec non_finite_node;
  non_finite_node.num_files = 2;
  non_finite_node.edit_nodes = [](std::uint64_t id, std::vector<Vector3>& nodes)
  {
    if (id == 5)
      nodes[0].x = std::numeric_limits<double>::infinity();
  };
  non_finite_node.edit_value = [](std::uint64_t /*id*/, double value)
  { return std::isfinite(value) ? value : 0.0; };
  ExpectRejected("non_finite_node_", non_finite_node);

  // A set cannot mix files with and without provenance.
  SetSpec untagged;
  untagged.num_files = 2;
  untagged.provenance_of = [](std::uint64_t file)
  {
    return file == 0 ? std::optional(std::make_pair(
                         static_cast<unsigned int>(FluxSolutionRole::FORWARD),
                         static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED)))
                     : std::nullopt;
  };
  ExpectRejected("untagged_", untagged);

  // An adjoint file cannot replace the flux moments of a forward problem.
  SetSpec adjoint;
  adjoint.num_files = 2;
  adjoint.role = static_cast<unsigned int>(FluxSolutionRole::ADJOINT);
  adjoint.orientation = static_cast<unsigned int>(AdjointFluxOrientation::SOLVER);
  const auto base = Write("adjoint_", adjoint);
  const auto phi_old = problem_->GetPhiOldLocal();
  EXPECT_THROW(LBSSolverIO::ReadFluxMoments(*problem_, base, false), std::runtime_error);
  EXPECT_EQ(problem_->GetPhiOldLocal(), phi_old);
}

TEST_F(CellwiseFileSetTest, RejectsInconsistentAngularSetsAndLeavesTheDestination)
{
  const std::vector<std::vector<double>> untouched = {{-1.0}};
  const auto ExpectRejected = [&](const std::string& name, const SetSpec& spec)
  {
    const auto base = Write(name, spec);
    auto psi = untouched;
    EXPECT_THROW(DiscreteOrdinatesProblemIO::ReadAngularFluxes(*problem_, base, psi),
                 std::runtime_error)
      << name;
    EXPECT_EQ(psi, untouched) << name;
  };

  // Adjoint files in both orientations, read without a requested orientation.
  SetSpec orientations;
  orientations.angular = true;
  orientations.num_files = 2;
  orientations.provenance_of = [](std::uint64_t file)
  {
    return std::optional(
      std::make_pair(static_cast<unsigned int>(FluxSolutionRole::ADJOINT),
                     static_cast<unsigned int>(file == 0 ? AdjointFluxOrientation::SOLVER
                                                         : AdjointFluxOrientation::PHYSICAL)));
  };
  ExpectRejected("psi_orientations_", orientations);

  // Forward and adjoint files.
  SetSpec roles;
  roles.angular = true;
  roles.num_files = 2;
  roles.provenance_of = [](std::uint64_t file)
  {
    return file == 0
             ? std::optional(
                 std::make_pair(static_cast<unsigned int>(FluxSolutionRole::FORWARD),
                                static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED)))
             : std::optional(
                 std::make_pair(static_cast<unsigned int>(FluxSolutionRole::ADJOINT),
                                static_cast<unsigned int>(AdjointFluxOrientation::PHYSICAL)));
  };
  ExpectRejected("psi_roles_", roles);

  // Directions without weights cannot be checked against the problem's quadrature.
  SetSpec no_weights;
  no_weights.angular = true;
  no_weights.num_files = 2;
  no_weights.omit_weights = true;
  ExpectRejected("psi_no_weights_", no_weights);
}
