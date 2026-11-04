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
#include "framework/utils/collective_error.h"
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
  // A failure on one rank is raised on every rank, so no rank goes on to build the problem alone.
  std::shared_ptr<MultiGroupXS> xs;
  ThrowIfAnyRankFailed(
    mpi_comm,
    CaptureError(
      [&]
      {
        {
          std::ofstream out(path);
          out << "NUM_GROUPS 2\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n0 1.0\n1 2.0\nSIGMA_T_END\n";
          out.close();
          if (out.fail())
            throw std::runtime_error("Failed to write " + path.string() + ".");
        }
        xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
      }),
    "Creating the test cross sections");
  return xs;
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
  /// Applied to the stored node locations of each written cell. The values stay those of the
  /// unedited nodes, so a node moved within the matching tolerance still reads exactly.
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
  /// Omit the quadrature directions and weights of each groupset of an angular set.
  bool omit_quadrature = false;
  /// The file_set_id attribute of each file, or nullopt to omit it; the default is one shared id.
  std::function<std::optional<std::uint64_t>(std::uint64_t)> file_set_id_of;
  /// The num_local_nodes attribute of every file, overriding the files' node counts.
  std::optional<std::uint64_t> recorded_nodes;
};

class CellwiseFileSetTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    // Under rank 0's working directory, unlike a node-local temporary directory, since rank 0
    // writes files that other ranks read. Rank 0 sends its full path, so every rank uses the same
    // directory even if the ranks were started in different working directories. A failure on
    // rank 0 is raised on every rank before the broadcast, so no rank waits in it or goes on to
    // build the problem.
    std::vector<std::string> path;
    RunWithCollectiveErrors(
      [&]
      {
        if (mpi_comm.rank() != 0)
          return;
        const auto nonce =
          static_cast<std::uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count());
        path.push_back(
          (std::filesystem::current_path() / ("opensn_cellwise_file_set_" + std::to_string(nonce)))
            .string());
        std::filesystem::create_directories(path.front());
      });
    BroadcastStrings(mpi_comm, path, 0);
    directory_ = path.front();
    problem_ = MakeProblem(directory_);
  }

  void TearDown() override
  {
    mpi_comm.barrier();
    if (mpi_comm.rank() == 0)
    {
      std::error_code error;
      std::filesystem::remove_all(directory_, error);
    }
  }

  /// Write a file set described by `spec` on rank 0 and return its base. Collective.
  std::string Write(const std::string& name, const SetSpec& spec)
  {
    const auto base = (directory_ / name).string();
    const auto cells = GatherCells(*problem_);
    // A failure on rank 0 is raised on every rank, and every file exists before any rank reads.
    RunWithCollectiveErrors(
      [&]
      {
        if (mpi_comm.rank() != 0)
          return;
        for (std::uint64_t file_id = 0; file_id < spec.num_files; ++file_id)
          WriteFile(spec, base, file_id, cells);
      });
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
    // The cells of this file, nodes reversed in every third cell, so node order differs. The
    // values are those of the nodes before `edit_nodes`.
    std::vector<CellNodes> file_cells;
    std::vector<std::vector<Vector3>> value_nodes;
    for (const auto& cell : cells)
    {
      const auto files = spec.files_of ? spec.files_of(cell.id)
                                       : std::vector<std::uint64_t>{(7 * cell.id) % spec.num_files};
      if (std::find(files.begin(), files.end(), file_id) == files.end())
        continue;
      auto& file_cell = file_cells.emplace_back(cell);
      if (cell.id % 3 == 0)
        std::reverse(file_cell.nodes.begin(), file_cell.nodes.end());
      value_nodes.push_back(file_cell.nodes);
      if (spec.edit_nodes)
        spec.edit_nodes(cell.id, file_cell.nodes);
    }
    if (spec.unknown_cell != 0 and file_id == 0)
    {
      file_cells.push_back({spec.unknown_cell, cells.front().nodes});
      value_nodes.push_back(cells.front().nodes);
    }

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
      for (std::size_t c = 0; c < file_cells.size(); ++c)
        for (const auto& node : value_nodes[c])
          for (std::size_t a = 0; a < num_unknowns; ++a)
            for (std::size_t g = 0; g < num_groups; ++g)
            {
              const auto value = Value(node, a, g, gs);
              values.push_back(spec.edit_value ? spec.edit_value(file_cells[c].id, value) : value);
            }
      return values;
    };

    // A failed fixture write throws, so that Write raises it on every rank before any rank reads.
    const auto require = [](const bool written)
    {
      if (not written)
        throw std::runtime_error("Failed to write a test file.");
    };
    const auto file_name = spec.single_file ? base + ".h5" : PartitionFileName(base, file_id);
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    require(file.Id() >= 0);
    const auto id = file.Id();
    const auto recorded =
      spec.recorded_files ? spec.recorded_files(file_id) : std::optional(spec.num_files);
    if (recorded.has_value())
    {
      require(H5CreateAttribute(id, "num_partitions", *recorded));
    }
    const auto file_set_id =
      spec.file_set_id_of ? spec.file_set_id_of(file_id) : std::optional<std::uint64_t>(7);
    if (file_set_id.has_value())
    {
      require(H5CreateAttribute(id, "file_set_id", *file_set_id));
    }
    const auto provenance = spec.provenance_of
                              ? spec.provenance_of(file_id)
                              : std::optional(std::make_pair(spec.role, spec.orientation));
    if (provenance.has_value())
    {
      require(H5CreateAttribute(id, "solution_role", provenance->first));
      require(H5CreateAttribute(id, "adjoint_orientation", provenance->second));
    }
    require(H5CreateGroup(id, "mesh"));
    require(H5CreateAttribute(id, "mesh/num_local_cells", std::uint64_t{cell_ids.size()}));
    require(H5CreateAttribute(
      id, "mesh/num_local_nodes", spec.recorded_nodes.value_or(std::uint64_t{nodes_x.size()})));
    require(H5WriteDataset1D(id, "mesh/cell_ids", cell_ids));
    require(H5WriteDataset1D(id, "mesh/num_cell_nodes", num_cell_nodes));
    require(H5WriteDataset1D(id, "mesh/nodes_x", nodes_x));
    require(H5WriteDataset1D(id, "mesh/nodes_y", nodes_y));
    require(H5WriteDataset1D(id, "mesh/nodes_z", nodes_z));

    if (not spec.angular)
    {
      require(H5CreateAttribute(id, "num_moments", problem_->GetNumMoments()));
      require(H5CreateAttribute(id, "num_groups", 2U));
      require(H5WriteDataset1D(id, "values", Values(problem_->GetNumMoments(), 2, 0)));
      return;
    }
    const auto& groupsets = problem_->GetGroupsets();
    require(H5CreateAttribute(id, "num_groupsets", std::uint64_t{groupsets.size()}));
    for (std::size_t gs = 0; gs < groupsets.size(); ++gs)
    {
      const auto& quadrature = groupsets[gs].quadrature;
      const auto group = "groupset_" + std::to_string(groupsets[gs].id);
      std::vector<double> omegas;
      for (const auto& omega : quadrature->GetOmegas())
        omegas.insert(omegas.end(), {omega.x, omega.y, omega.z});
      require(H5CreateGroup(id, group));
      require(H5CreateAttribute(
        id, group + "/num_directions", std::uint64_t{quadrature->GetNumAngles()}));
      require(H5CreateAttribute(id, group + "/num_groups", 1U));
      if (not spec.omit_quadrature)
      {
        require(H5WriteDataset1D(id, group + "/omegas", omegas));
      }
      if (not spec.omit_weights and not spec.omit_quadrature)
      {
        require(H5WriteDataset1D(id, group + "/weights", quadrature->GetWeights()));
      }
      require(H5WriteDataset1D(id, group + "/values", Values(quadrature->GetNumAngles(), 1, gs)));
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
  // The write is collective, so every file exists before any rank reads.
  LBSSolverIO::WriteFluxMoments(*problem_, base);
  std::vector<double> read;
  LBSSolverIO::ReadFluxMoments(*problem_, base, false, read);
  EXPECT_EQ(read, phi);
}

TEST_F(CellwiseFileSetTest, RaisesAWriteFailureOnEveryRank)
{
  // A directory where rank P - 1's file belongs makes only that rank's write fail.
  const auto base = (directory_ / "unwritable_").string();
  RunWithCollectiveErrors(
    [&]
    {
      if (mpi_comm.rank() == 0)
        std::filesystem::create_directories(PartitionFileName(base, mpi_comm.size() - 1));
    });
  EXPECT_THROW(LBSSolverIO::WriteFluxMoments(*problem_, base), std::runtime_error);
  // The test problem does not save angular fluxes, so a source of the right size is supplied.
  std::vector<std::vector<double>> psi;
  for (const auto& groupset : problem_->GetGroupsets())
    psi.emplace_back(problem_->GetSpatialDiscretization().GetNumLocalDOFs(groupset.psi_uk_man_),
                     0.0);
  EXPECT_THROW(DiscreteOrdinatesProblemIO::WriteAngularFluxes(*problem_, base, psi),
               std::runtime_error);
  // A source of the wrong size is rejected on every rank before any file is written.
  std::vector<std::vector<double>> wrong(psi.size());
  EXPECT_THROW(DiscreteOrdinatesProblemIO::WriteAngularFluxes(
                 *problem_, (directory_ / "wrong_size_").string(), wrong),
               std::runtime_error);
  // Every rank left the writes, so a following collective completes.
  int count = 1;
  mpi_comm.all_reduce(count, mpi::op::sum<int>());
  EXPECT_EQ(count, mpi_comm.size());
}

TEST_F(CellwiseFileSetTest, RejectedWriteLeavesTheSetUntouched)
{
  auto& phi = problem_->GetPhiNewLocal();
  for (std::size_t k = 0; k < phi.size(); ++k)
    phi[k] = 0.25 * static_cast<double>(k % 97);
  const auto base = (directory_ / "kept_").string();
  LBSSolverIO::WriteFluxMoments(*problem_, base);

  // Rank P - 1 offers a source of the wrong size and the other ranks valid sources with other
  // values. Every rank checks its source before any rank replaces its file, so the call fails on
  // every rank and the set still holds the first write.
  std::vector<double> other(phi.size(), 1.0);
  if (mpi_comm.rank() == mpi_comm.size() - 1)
    other.push_back(1.0);
  EXPECT_THROW(LBSSolverIO::WriteFluxMoments(*problem_, base, other), std::runtime_error);
  std::vector<double> read;
  LBSSolverIO::ReadFluxMoments(*problem_, base, false, read);
  EXPECT_EQ(read, phi);

  // The angular writer checks its source separately, before any rank replaces its file.
  std::vector<std::vector<double>> psi;
  for (const auto& groupset : problem_->GetGroupsets())
  {
    auto& values =
      psi.emplace_back(problem_->GetSpatialDiscretization().GetNumLocalDOFs(groupset.psi_uk_man_));
    for (std::size_t k = 0; k < values.size(); ++k)
      values[k] = 0.5 * static_cast<double>(k % 89);
  }
  const auto psi_base = (directory_ / "kept_psi_").string();
  DiscreteOrdinatesProblemIO::WriteAngularFluxes(*problem_, psi_base, psi);

  auto other_psi = psi;
  for (auto& values : other_psi)
    std::fill(values.begin(), values.end(), 1.0);
  if (mpi_comm.rank() == mpi_comm.size() - 1)
    other_psi.back().push_back(1.0);
  EXPECT_THROW(DiscreteOrdinatesProblemIO::WriteAngularFluxes(*problem_, psi_base, other_psi),
               std::runtime_error);
  std::vector<std::vector<double>> read_psi;
  DiscreteOrdinatesProblemIO::ReadAngularFluxes(*problem_, psi_base, read_psi);
  EXPECT_EQ(read_psi, psi);
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
  // The set is written outside the expectation, so a failure to write it cannot satisfy it.
  const auto more_base = Write("unrecorded_more_", spec);
  EXPECT_THROW(LBSSolverIO::ReadFluxMoments(*problem_, more_base, false, phi), std::runtime_error);
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
  const auto of_3_base = Write("single_of_3", spec);
  EXPECT_THROW(LBSSolverIO::ReadFluxMoments(*problem_, of_3_base, true, phi), std::runtime_error);
}

TEST_F(CellwiseFileSetTest, MatchesNodesWithinTheTolerance)
{
  // One node of a cell stored a few rounding errors away from the mesh node matches it through the
  // tolerance, while the other nodes of the cell match exactly.
  SetSpec spec;
  spec.num_files = 2;
  spec.edit_nodes = [](std::uint64_t id, std::vector<Vector3>& nodes)
  {
    if (id == 5)
      nodes[0].x += 8.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, nodes[0].x);
  };
  std::vector<double> phi;
  LBSSolverIO::ReadFluxMoments(*problem_, Write("near_node_", spec), false, phi);
  ExpectFluxMoments(phi);
}

TEST_F(CellwiseFileSetTest, RejectsAnImpossibleFileCountUnread)
{
  // A recorded count above INT_MAX is rejected on rank 0 before any file is validated or read.
  SetSpec spec;
  spec.num_files = 1;
  spec.recorded_files = [](std::uint64_t)
  { return std::optional(std::numeric_limits<std::uint64_t>::max()); };
  const auto base = Write("impossible_", spec);
  int num_validated = 0;
  EXPECT_THROW(ReadCellwiseFileSet(
                 base,
                 false,
                 {{"values", 2ULL * problem_->GetNumMoments()}},
                 *problem_->GetGrid(),
                 problem_->GetSpatialDiscretization(),
                 [&](hid_t, const std::string&) { ++num_validated; },
                 [](std::size_t, const Cell&, std::size_t, const double*) {},
                 mpi_comm,
                 "Reading test flux moments"),
               std::runtime_error);
  EXPECT_EQ(num_validated, 0);

  // A zero send budget is rejected before any file is opened.
  EXPECT_THROW(ReadCellwiseFileSet(
                 base,
                 false,
                 {{"values", 2ULL * problem_->GetNumMoments()}},
                 *problem_->GetGrid(),
                 problem_->GetSpatialDiscretization(),
                 [&](hid_t, const std::string&) { ++num_validated; },
                 [](std::size_t, const Cell&, std::size_t, const double*) {},
                 mpi_comm,
                 "Reading test flux moments",
                 0),
               std::invalid_argument);
  EXPECT_EQ(num_validated, 0);
}

TEST_F(CellwiseFileSetTest, ReturnsTheRecordedOrientationOnEveryRank)
{
  // With one file, the ranks without a file take the set's provenance from the others.
  SetSpec spec;
  spec.num_files = 1;
  spec.role = static_cast<unsigned int>(FluxSolutionRole::ADJOINT);
  spec.orientation = static_cast<unsigned int>(AdjointFluxOrientation::PHYSICAL);
  std::vector<double> phi;
  EXPECT_EQ(LBSSolverIO::ReadFluxMoments(*problem_, Write("physical_", spec), false, phi),
            AdjointFluxOrientation::PHYSICAL);
  ExpectFluxMoments(phi);

  // A set without provenance records no orientation.
  spec.provenance_of = [](std::uint64_t)
  { return std::optional<std::pair<unsigned int, unsigned int>>(); };
  EXPECT_EQ(LBSSolverIO::ReadFluxMoments(*problem_, Write("untagged_set_", spec), false, phi),
            AdjointFluxOrientation::UNSPECIFIED);
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

  // An infinite coordinate would make the node-matching tolerance infinite. The values stay those
  // of the unedited nodes, so only the coordinate is invalid.
  SetSpec non_finite_node;
  non_finite_node.num_files = 2;
  non_finite_node.edit_nodes = [](std::uint64_t id, std::vector<Vector3>& nodes)
  {
    if (id == 5)
      nodes[0].x = std::numeric_limits<double>::infinity();
  };
  ExpectRejected("non_finite_node_", non_finite_node);

  // A finite node so far away that its squared distance to any mesh node overflows to infinity.
  SetSpec far_node;
  far_node.num_files = 2;
  far_node.edit_nodes = [](std::uint64_t id, std::vector<Vector3>& nodes)
  {
    if (id == 5)
      nodes[0].x = 1.0e200;
  };
  ExpectRejected("far_node_", far_node);

  // A recorded file count that no MPI run can have written. (RejectsAnImpossibleFileCountUnread
  // checks that it is rejected before any file is read.)
  SetSpec huge_count;
  huge_count.num_files = 2;
  huge_count.recorded_files = [](std::uint64_t /*file*/)
  { return std::optional(std::numeric_limits<std::uint64_t>::max()); };
  ExpectRejected("huge_count_", huge_count);

  // Files of two writes, as after an interrupted overwrite, and files that record their number but
  // not their write.
  SetSpec mixed_write;
  mixed_write.num_files = 2;
  mixed_write.file_set_id_of = [](std::uint64_t file)
  { return std::optional<std::uint64_t>(file == 0 ? 7 : 8); };
  ExpectRejected("mixed_write_", mixed_write);
  SetSpec no_set_id;
  no_set_id.num_files = 2;
  no_set_id.file_set_id_of = [](std::uint64_t) { return std::optional<std::uint64_t>(); };
  ExpectRejected("no_set_id_", no_set_id);

  // A file that records more nodes than the whole mesh has.
  SetSpec too_many_nodes;
  too_many_nodes.num_files = 2;
  too_many_nodes.recorded_nodes = std::uint64_t{1} << 40;
  ExpectRejected("too_many_nodes_", too_many_nodes);

  // Two nodes of a cell moved within the tolerance of the same mesh node: neither matches exactly,
  // and they share their only candidate, so the match is ambiguous.
  SetSpec ambiguous;
  ambiguous.num_files = 2;
  ambiguous.edit_nodes = [](std::uint64_t id, std::vector<Vector3>& nodes)
  {
    if (id != 5)
      return;
    const auto shift = 8.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, nodes[0].x);
    nodes[1] = nodes[0];
    nodes[0].x += shift;
    nodes[1].y += shift;
  };
  ExpectRejected("ambiguous_", ambiguous);

  // File 0 does not record the number of files, but a later file does. With one rank the set
  // has only file 0.
  if (mpi_comm.size() > 1)
  {
    SetSpec mixed_reverse;
    mixed_reverse.num_files = static_cast<std::uint64_t>(mpi_comm.size());
    mixed_reverse.recorded_files = [](std::uint64_t file)
    {
      return file == 0 ? std::nullopt
                       : std::optional<std::uint64_t>(static_cast<std::uint64_t>(mpi_comm.size()));
    };
    ExpectRejected("mixed_reverse_", mixed_reverse);
  }

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

  // Nor files that record an unspecified role, as a file written from a supplied vector does, and
  // files without provenance.
  SetSpec unspecified_mixed;
  unspecified_mixed.num_files = 2;
  unspecified_mixed.provenance_of = [](std::uint64_t file)
  {
    return file == 0 ? std::optional(std::make_pair(
                         static_cast<unsigned int>(FluxSolutionRole::UNSPECIFIED),
                         static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED)))
                     : std::nullopt;
  };
  ExpectRejected("unspecified_mixed_", unspecified_mixed);

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

  // Files that record an unspecified role and files without provenance.
  SetSpec unspecified_mixed;
  unspecified_mixed.angular = true;
  unspecified_mixed.num_files = 2;
  unspecified_mixed.provenance_of = [](std::uint64_t file)
  {
    return file == 0 ? std::optional(std::make_pair(
                         static_cast<unsigned int>(FluxSolutionRole::UNSPECIFIED),
                         static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED)))
                     : std::nullopt;
  };
  ExpectRejected("psi_unspecified_mixed_", unspecified_mixed);

  // Directions without weights cannot be checked against the problem's quadrature.
  SetSpec no_weights;
  no_weights.angular = true;
  no_weights.num_files = 2;
  no_weights.omit_weights = true;
  ExpectRejected("psi_no_weights_", no_weights);

  // A file that records the number of files in its set is never written without its quadrature.
  SetSpec no_quadrature;
  no_quadrature.angular = true;
  no_quadrature.num_files = 2;
  no_quadrature.omit_quadrature = true;
  ExpectRejected("psi_no_quadrature_", no_quadrature);

  // Files of two writes.
  SetSpec mixed_write;
  mixed_write.angular = true;
  mixed_write.num_files = 2;
  mixed_write.file_set_id_of = [](std::uint64_t file)
  { return std::optional<std::uint64_t>(file == 0 ? 7 : 8); };
  ExpectRejected("psi_mixed_write_", mixed_write);
}
