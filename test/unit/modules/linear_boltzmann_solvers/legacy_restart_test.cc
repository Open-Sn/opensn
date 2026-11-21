// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include "framework/utils/hdf_utils.h"
#include <gtest/gtest.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

constexpr unsigned int kGroups = 3;

std::shared_ptr<MultiGroupXS>
MakeThreeGroupXS()
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_legacy_restart_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << "NUM_GROUPS 3\nNUM_MOMENTS 1\nSIGMA_T_BEGIN\n0 1.0\n1 2.0\n2 3.0\nSIGMA_T_END\n";
  }
  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
  std::filesystem::remove(path);
  return xs;
}

/// A 2D problem with P1 moments and two groupsets of unequal sizes, (0) and (1, 2), so the solver
/// stores flux moments by groupset rather than in the [node][moment][group] order of legacy files.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem()
{
  auto grid = BuildSquareMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  const auto quadrature =
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature2DXY>(2, 4, 1));
  ParameterBlock groupsets("groupsets");
  const std::vector<std::vector<int>> ranges = {{0, 0}, {1, 2}};
  for (std::size_t gs = 0; gs < ranges.size(); ++gs)
  {
    ParameterBlock groupset(std::to_string(gs));
    groupset.AddParameter(ParameterBlock("groups_from_to", ranges[gs]));
    groupset.AddParameter("angular_quadrature", quadrature);
    groupsets.AddParameter(groupset);
  }
  groupsets.ChangeToArray();

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeThreeGroupXS());
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", static_cast<int>(kGroups));
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  return DiscreteOrdinatesProblem::Create(params);
}

/// Distinct values in [node][moment][group] order, the order of a legacy file.
std::vector<double>
MakeLegacyValues(std::size_t size, double offset)
{
  std::vector<double> values(size);
  for (std::size_t k = 0; k < size; ++k)
    values[k] = offset + static_cast<double>(k);
  return values;
}

/// This rank's file of a restart named by a stem that every rank shares.
std::filesystem::path
RankPath(const std::string& name)
{
  std::int64_t stamp = std::chrono::steady_clock::now().time_since_epoch().count();
  mpi_comm.broadcast(stamp, 0);
  return std::filesystem::temp_directory_path() / (name + "_" + std::to_string(stamp) + "_" +
                                                   std::to_string(mpi_comm.rank()) + ".restart.h5");
}

/// Writes a version-1 restart file: raw local vectors in [node][moment][group] order.
void
WriteLegacyRestart(const std::filesystem::path& path,
                   const std::vector<double>& phi_old,
                   const std::vector<double>* phi_new)
{
  const H5FileHandle file(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  ASSERT_GE(file.Id(), 0);
  ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), "phi_old", phi_old));
  if (phi_new != nullptr)
  {
    ASSERT_TRUE(H5WriteDataset1D<double>(file.Id(), "phi_new", *phi_new));
  }
}

/// Expects `stored`, in the problem's layout, to hold `legacy` given in [node][moment][group]
/// order.
void
ExpectLegacyOrder(const DiscreteOrdinatesProblem& problem,
                  std::span<const double> stored,
                  const std::vector<double>& legacy)
{
  const auto& layout = problem.GetFluxMomentLayout();
  ASSERT_EQ(stored.size(), legacy.size());
  for (std::uint64_t node = 0; node < layout.GetNumLocalNodes(); ++node)
    for (unsigned int m = 0; m < layout.GetNumMoments(); ++m)
      for (unsigned int g = 0; g < kGroups; ++g)
        EXPECT_EQ(stored[layout.MapDOF(node, m, g)],
                  legacy[(node * layout.GetNumMoments() + m) * kGroups + g]);
}

} // namespace

TEST(LegacyRestartTest, RawFluxIsReorderedIntoTheGroupsetLayout)
{
  auto problem = MakeProblem();
  const auto& layout = problem->GetFluxMomentLayout();
  ASSERT_FALSE(layout.IsNodeMomentGroupOrder());
  ASSERT_GT(layout.GetNumMoments(), 1U);

  const auto legacy_old = MakeLegacyValues(layout.GetNumLocalDOFs(), 1.0);
  const auto legacy_new = MakeLegacyValues(layout.GetNumLocalDOFs(), -1.0e3);
  const auto path = RankPath("opensn_legacy_restart");
  WriteLegacyRestart(path, legacy_old, &legacy_new);
  mpi_comm.barrier();

  EXPECT_TRUE(problem->ReadRestartData({}, path));
  ExpectLegacyOrder(*problem, std::as_const(*problem).GetPhiOldLocal().Data(), legacy_old);
  ExpectLegacyOrder(*problem, std::as_const(*problem).GetPhiNewLocal().Data(), legacy_new);

  mpi_comm.barrier();
  std::filesystem::remove(path);
}

TEST(LegacyRestartTest, MissingNewFluxIsTheOldFlux)
{
  auto problem = MakeProblem();
  const auto& layout = problem->GetFluxMomentLayout();
  const auto legacy_old = MakeLegacyValues(layout.GetNumLocalDOFs(), 0.5);
  const auto path = RankPath("opensn_legacy_restart_old_only");
  WriteLegacyRestart(path, legacy_old, nullptr);
  mpi_comm.barrier();

  EXPECT_TRUE(problem->ReadRestartData({}, path));
  ExpectLegacyOrder(*problem, std::as_const(*problem).GetPhiOldLocal().Data(), legacy_old);
  ExpectLegacyOrder(*problem, std::as_const(*problem).GetPhiNewLocal().Data(), legacy_old);

  mpi_comm.barrier();
  std::filesystem::remove(path);
}
