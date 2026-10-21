// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "test/unit/common/mesh_builders.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/uncollided_problem/uncollided_problem.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

using namespace opensn;

namespace
{

// Two groups with P1 scattering, up- and downscattering, and fission into both groups. All values
// are dyadic, so the source moments made from dyadic flux moments are exact in double precision.
// kTransfer[ell][g][gp] is the transfer from group gp to group g, kProduction[g][gp] = chi[g] *
// nu * sigma_f[gp] with chi = (0.75, 0.25), nu = 2.5 and sigma_f = (0.0625, 0.125).
constexpr double kTransfer[2][2][2] = {{{0.5, 0.125}, {0.25, 1.0}}, {{0.375, 0.0}, {0.0625, 0.5}}};
constexpr double kProduction[2][2] = {{0.1171875, 0.234375}, {0.0390625, 0.078125}};

std::shared_ptr<MultiGroupXS>
LoadXS(const std::string& contents)
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_source_moments_from_phi_r" + std::to_string(mpi_comm.rank()) + "_" +
                     std::to_string(now) + ".xs");
  {
    std::ofstream out(path);
    out << contents;
  }
  auto xs = std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(path.string()));
  std::filesystem::remove(path);
  return xs;
}

std::shared_ptr<MultiGroupXS>
MakeTwoGroupXS()
{
  return LoadXS("NUM_GROUPS 2\nNUM_MOMENTS 2\n"
                "SIGMA_T_BEGIN\n0 1.0\n1 2.0\nSIGMA_T_END\n"
                "SIGMA_F_BEGIN\n0 0.0625\n1 0.125\nSIGMA_F_END\n"
                "NU_BEGIN\n0 2.5\n1 2.5\nNU_END\n"
                "CHI_BEGIN\n0 0.75\n1 0.25\nCHI_END\n"
                "TRANSFER_MOMENTS_BEGIN\n"
                "M_GFROM_GTO_VAL 0 0 0 0.5\n"
                "M_GFROM_GTO_VAL 0 0 1 0.25\n"
                "M_GFROM_GTO_VAL 0 1 0 0.125\n"
                "M_GFROM_GTO_VAL 0 1 1 1.0\n"
                "M_GFROM_GTO_VAL 1 0 0 0.375\n"
                "M_GFROM_GTO_VAL 1 0 1 0.0625\n"
                "M_GFROM_GTO_VAL 1 1 1 0.5\n"
                "TRANSFER_MOMENTS_END\n");
}

/// A 2D XY problem with P1 scattering and one groupset per group, so that the source moments
/// include across-groupset scattering and fission.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem()
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
  xs_entry.AddParameter("xs", MakeTwoGroupXS());
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

/// Dyadic values offset + step * (k % period), k the local index.
std::vector<double>
MakeVector(std::size_t size, double offset, double step, std::size_t period)
{
  std::vector<double> values(size);
  for (std::size_t k = 0; k < size; ++k)
    values[k] = offset + step * static_cast<double>(k % period);
  return values;
}

/// The source moments of the cross sections above: at every node and moment of degree ell,
/// q_g = sum_gp kTransfer[ell][g][gp] phi_gp, plus sum_gp kProduction[g][gp] phi_gp for ell = 0.
std::vector<double>
ReferenceSourceMoments(const LBSProblem& problem, const std::vector<double>& phi)
{
  const auto& sdm = problem.GetSpatialDiscretization();
  const auto& uk_man = problem.GetUnknownManager();
  const auto& moment_map = problem.GetGroupset(0).quadrature->GetMomentToHarmonicsIndexMap();
  std::vector<double> q(phi.size(), 0.0);
  for (const auto& cell : problem.GetGrid()->GetLocalCells())
    for (unsigned int i = 0; i < sdm.GetCellNumNodes(*cell); ++i)
      for (unsigned int m = 0; m < problem.GetNumMoments(); ++m)
      {
        const auto ell = moment_map[m].ell;
        for (unsigned int g = 0; g < 2; ++g)
        {
          double value = 0.0;
          for (unsigned int gp = 0; gp < 2; ++gp)
          {
            const double phi_gp = phi[sdm.MapDOFLocal(*cell, i, uk_man, m, gp)];
            value += kTransfer[ell][g][gp] * phi_gp;
            if (ell == 0)
              value += kProduction[g][gp] * phi_gp;
          }
          q[sdm.MapDOFLocal(*cell, i, uk_man, m, g)] = value;
        }
      }
  return q;
}

} // namespace

// The form without an argument makes the source moments from phi_new: with phi_old different from
// phi_new, it gives bit for bit the result of the overload given phi_new, equal to the reference
// source moments of phi_new and different from those of phi_old.
TEST(MakeSourceMomentsFromPhiTest, NoArgumentFormUsesPhiNew)
{
  auto problem = MakeProblem();
  const auto num_local_dofs = problem->GetPhiNewLocal().size();
  ASSERT_EQ(num_local_dofs,
            problem->GetSpatialDiscretization().GetNumLocalDOFs(problem->GetUnknownManager()));
  ASSERT_EQ(problem->GetNumMoments(), 3U);

  const auto phi_new = MakeVector(num_local_dofs, 1.0, 0.125, 11);
  const auto phi_old = MakeVector(num_local_dofs, 4.0, -0.25, 5);
  problem->SetPhiNewFrom(phi_new);
  problem->SetPhiOldFrom(phi_old);

  const auto from_phi_new = problem->MakeSourceMomentsFromPhi();
  const auto expected = ReferenceSourceMoments(*problem, phi_new);
  EXPECT_EQ(from_phi_new, problem->MakeSourceMomentsFromPhi(phi_new));
  EXPECT_EQ(from_phi_new, expected);
  EXPECT_NE(from_phi_new, ReferenceSourceMoments(*problem, phi_old));
}

// The overload makes the source moments from the given flux moments and leaves phi_new, phi_old and
// the external source moments unchanged.
TEST(MakeSourceMomentsFromPhiTest, OverloadLeavesProblemStateUnchanged)
{
  auto problem = MakeProblem();
  const auto num_local_dofs = problem->GetPhiNewLocal().size();

  const auto phi_new = MakeVector(num_local_dofs, 1.0, 0.125, 11);
  const auto phi_old = MakeVector(num_local_dofs, 4.0, -0.25, 5);
  const auto ext_src = MakeVector(num_local_dofs, -1.0, -0.5, 3);
  const auto phi = MakeVector(num_local_dofs, 0.5, 0.0625, 7);
  problem->SetPhiNewFrom(phi_new);
  problem->SetPhiOldFrom(phi_old);
  problem->SetExtSrcMomentsFrom(ext_src);

  const auto source_moments = problem->MakeSourceMomentsFromPhi(phi);
  EXPECT_EQ(source_moments, ReferenceSourceMoments(*problem, phi));
  EXPECT_NE(source_moments, ReferenceSourceMoments(*problem, phi_new));
  EXPECT_EQ(problem->GetPhiNewLocal(), phi_new);
  EXPECT_EQ(problem->GetPhiOldLocal(), phi_old);
  EXPECT_EQ(problem->GetExtSrcMomentsLocal(), ext_src);
}

// Flux moments whose size is not the problem's number of local DOFs are refused before any state
// changes. The short vector has the capacity of a correct one, so that a missing check reads
// allocated memory, and it is tried only after the long one has been refused.
TEST(MakeSourceMomentsFromPhiTest, WrongSizeThrows)
{
  auto problem = MakeProblem();
  const auto num_local_dofs = problem->GetPhiNewLocal().size();
  ASSERT_GT(num_local_dofs, 0U);

  const auto phi_new = MakeVector(num_local_dofs, 1.0, 0.125, 11);
  const auto phi_old = MakeVector(num_local_dofs, 4.0, -0.25, 5);
  problem->SetPhiNewFrom(phi_new);
  problem->SetPhiOldFrom(phi_old);

  const auto too_long = MakeVector(num_local_dofs + 1, 0.5, 0.0625, 7);
  std::vector<double> too_short;
  too_short.reserve(num_local_dofs);
  too_short.assign(num_local_dofs - 1, 0.5);
  ASSERT_THROW(problem->MakeSourceMomentsFromPhi(too_long), std::invalid_argument);
  EXPECT_THROW(problem->MakeSourceMomentsFromPhi(too_short), std::invalid_argument);
  EXPECT_EQ(problem->GetPhiNewLocal(), phi_new);
  EXPECT_EQ(problem->GetPhiOldLocal(), phi_old);
}

// An UncollidedProblem has no source function, so making the source moments throws a
// std::logic_error (not a subclass such as std::invalid_argument) whose message names the problem
// and the operation, and the problem's flux moments and external source moments are left
// unchanged. The problem has its own name, so that the message is checked for the problem's name
// and not for its type.
TEST(MakeSourceMomentsFromPhiTest, UncollidedProblemThrowsWithoutStateChange)
{
  if (mpi_comm.size() != 1)
    GTEST_SKIP() << "Uncollided flux generation runs with exactly one MPI rank.";

  auto grid = BuildSquareMesh(1.0, 4, 0.0);
  grid->SetUniformBlockID(0);

  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 0}));
  ParameterBlock groupsets("groupsets");
  groupsets.AddParameter(groupset);
  groupsets.ChangeToArray();

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter(
    "xs", std::make_shared<MultiGroupXS>(MultiGroupXS::CreateSimpleOneGroup(1.0, 0.5)));
  ParameterBlock xs_map("xs_map");
  xs_map.AddParameter(xs_entry);
  xs_map.ChangeToArray();

  ParameterBlock point_source;
  point_source.AddParameter(ParameterBlock("location", std::vector<double>{0.3, 0.4, 0.0}));
  point_source.AddParameter(ParameterBlock("strength", std::vector<double>{1.0}));
  ParameterBlock near_source;
  near_source.AddParameter("xmin", -0.01);
  near_source.AddParameter("xmax", 1.01);
  near_source.AddParameter("ymin", -0.01);
  near_source.AddParameter("ymax", 1.01);
  near_source.AddParameter("infz", true);

  const std::string problem_name = "uncollided_point_source";
  ParameterBlock params;
  params.AddParameter("name", problem_name);
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 1);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter(ParameterBlock(
    "point_sources", std::vector<std::shared_ptr<PointSource>>{PointSource::Create(point_source)}));
  params.AddParameter(ParameterBlock(
    "near_source",
    std::vector<std::shared_ptr<LogicalVolume>>{RPPLogicalVolume::Create(near_source)}));
  auto problem = UncollidedProblem::Create(params);
  ASSERT_EQ(problem->GetName(), problem_name);

  const auto num_local_dofs = problem->GetPhiNewLocal().size();
  const auto phi_new = MakeVector(num_local_dofs, 1.0, 0.125, 11);
  const auto phi_old = MakeVector(num_local_dofs, 4.0, -0.25, 5);
  const auto ext_src = MakeVector(num_local_dofs, -1.0, -0.5, 3);
  problem->SetPhiNewFrom(phi_new);
  problem->SetPhiOldFrom(phi_old);
  problem->SetExtSrcMomentsFrom(ext_src);

  const auto phi = MakeVector(num_local_dofs, 0.5, 0.0625, 7);
  EXPECT_THROW(
    {
      try
      {
        problem->MakeSourceMomentsFromPhi(phi);
      }
      catch (const std::logic_error& error)
      {
        const std::string message = error.what();
        EXPECT_TRUE(typeid(error) == typeid(std::logic_error)) << typeid(error).name();
        EXPECT_NE(message.find(problem_name), std::string::npos) << message;
        EXPECT_NE(message.find("MakeSourceMomentsFromPhi"), std::string::npos) << message;
        throw;
      }
    },
    std::logic_error);
  EXPECT_EQ(problem->GetPhiNewLocal(), phi_new);
  EXPECT_EQ(problem->GetPhiOldLocal(), phi_old);
  EXPECT_EQ(problem->GetExtSrcMomentsLocal(), ext_src);
}
