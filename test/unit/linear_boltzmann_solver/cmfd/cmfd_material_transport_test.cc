// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_acceleration.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/pi_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

namespace
{
std::shared_ptr<DiscreteOrdinatesProblem>
MakeMixedSlab()
{
  auto grid = BuildLineMesh(4.0, 8, 0.0);
  grid->SetUniformBlockID(0);
  for (unsigned int c = 1; c < 8; c += 2)
  {
    ParameterBlock volume;
    volume.AddParameter("infx", true);
    volume.AddParameter("infy", true);
    volume.AddParameter("zmin", 0.5 * c);
    volume.AddParameter("zmax", 0.5 * (c + 1));
    grid->SetBlockIDFromLogicalVolume(*RPPLogicalVolume::Create(volume), 1, true);
  }
  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 1}));
  groupset.AddParameter(
    "angular_quadrature",
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLProductQuadrature1DSlab>(4, 0)));
  groupset.AddParameter("inner_linear_method", std::string("classic_richardson"));
  groupset.AddParameter("l_max_its", 1);
  groupset.AddParameter("l_abs_tol", 1.0e-14);
  ParameterBlock groupsets("groupsets");
  groupsets.ChangeToArray();
  groupsets.AddParameter(groupset);
  ParameterBlock xs_map("xs_map");
  xs_map.ChangeToArray();
  for (int material = 0; material < 2; ++material)
  {
    const std::string name = material == 0 ? "cmfd_material_a_2g.xs" : "cmfd_material_b_2g.xs";
    auto xs = std::make_shared<MultiGroupXS>(
      MultiGroupXS::LoadFromOpenSn(std::string(OPENSN_TEST_ROOT) + "/assets/xs/" + name));
    ParameterBlock entry(std::to_string(material));
    entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{material}));
    entry.AddParameter("xs", xs);
    xs_map.AddParameter(entry);
  }
  ParameterBlock boundary("0");
  boundary.AddParameter("name", std::string("zmin"));
  boundary.AddParameter("type", std::string("reflecting"));
  ParameterBlock boundaries("boundary_conditions");
  boundaries.ChangeToArray();
  boundaries.AddParameter(boundary);
  ParameterBlock options("options");
  options.AddParameter("verbose_inner_iterations", false);
  options.AddParameter("verbose_outer_iterations", false);
  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 2);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter(boundaries);
  params.AddParameter(options);
  return DiscreteOrdinatesProblem::Create(params);
}
} // namespace

TEST(CMFDMaterialTransport, MixedAggregatesMatchIndependentVacuumDGReference)
{
  // The reference is a separately assembled generalized eigensystem, not OpenSn output.
  // Its generator also verifies neutron balance and the algebraic eigenpair residual.
  std::ifstream input(std::string(OPENSN_TEST_ROOT) + "/assets/xs/cmfd_material_dg_reference.txt");
  ASSERT_TRUE(input.good());
  std::string header;
  std::getline(input, header);
  double expected_k = 0.0;
  std::array<double, 2> expected_currents{};
  std::array<double, 32> expected_phi{};
  input >> expected_k;
  for (auto& value : expected_currents)
    input >> value;
  for (auto& value : expected_phi)
    input >> value;
  ASSERT_FALSE(input.fail());
  ASSERT_GT(expected_k, 0.0);

  // The CMFD fixed point is the transport fixed point, so OpenSn must reproduce the PWLD
  // eigenpair that the generator assembles independently. The differences that remain are
  // iteration error (k_tol 1e-12, fission_source_tol 1e-10 and the converged coarse solves) and
  // round-off, far below the tolerances used here.
  for (const int energy_bin : {1, 2})
  {
    SCOPED_TRACE(::testing::Message() << "energy_bin=" << energy_bin);
    auto problem = MakeMixedSlab();
    ParameterBlock cmfd_params;
    cmfd_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
    cmfd_params.AddParameter("name", std::string("mixed_dg_"));
    ParameterBlock coarse_mesh("coarse_mesh");
    coarse_mesh.AddParameter("type", std::string("local_aggregation"));
    coarse_mesh.AddParameter("size", 2);
    cmfd_params.AddParameter(coarse_mesh);
    cmfd_params.AddParameter("coarse_groups", 2 / energy_bin);
    cmfd_params.AddParameter("petsc_options",
                             std::string("-mixed_dg_ksp_type gmres -mixed_dg_pc_type jacobi"));
    auto cmfd = CMFDAcceleration::Create(cmfd_params);
    ParameterBlock solver_params;
    solver_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
    solver_params.AddParameter("acceleration",
                               std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>(cmfd));
    solver_params.AddParameter("max_iters", 500);
    solver_params.AddParameter("k_tol", 1.0e-12);
    solver_params.AddParameter("fission_source_tol", 1.0e-10);
    auto solver = PowerIterationKEigenSolver::Create(solver_params);
    solver->Initialize();
    solver->Execute();
    EXPECT_TRUE(solver->IsConverged());
    EXPECT_GT(cmfd->GetNumAppliedCorrections(), 0);
    EXPECT_NEAR(solver->GetEigenvalue(), expected_k, 1.0e-9);
    const auto& phi = problem->GetPhiNewLocal();
    const double production = ComputeFissionProduction(*problem, phi);
    ASSERT_GT(production, 0.0);
    const auto vacuum_boundary = problem->FindBoundaryID("zmax");
    ASSERT_TRUE(vacuum_boundary.has_value());
    const auto leakage = ComputeLeakage(*problem, 0, vacuum_boundary.value_or(0));
    ASSERT_EQ(leakage.size(), 2);
    for (std::size_t g = 0; g < leakage.size(); ++g)
      EXPECT_NEAR(leakage[g] / production, expected_currents[g], 1.0e-8);
    const auto& sdm = problem->GetSpatialDiscretization();
    const auto& unknowns = problem->GetUnknownManager();
    for (const auto& cell : problem->GetGrid()->GetLocalCells())
    {
      const auto c = static_cast<std::size_t>(cell->centroid.z / 0.5);
      const auto& nodes = sdm.GetCellNodeLocations(*cell);
      ASSERT_LT(c, 8);
      ASSERT_EQ(nodes.size(), 2);
      for (unsigned int i = 0; i < nodes.size(); ++i)
      {
        const auto n = nodes[i].z < cell->centroid.z ? 0 : 1;
        for (unsigned int g = 0; g < 2; ++g)
          EXPECT_NEAR(phi[sdm.MapDOFLocal(*cell, i, unknowns, 0, g)] / production,
                      expected_phi[(c * 2 + n) * 2 + g],
                      1.0e-7);
      }
    }
  }
}
