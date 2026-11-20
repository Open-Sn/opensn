// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/mesh/mesh_generator/extruder_mesh_generator.h"
#include "framework/mesh/mesh_generator/from_file_mesh_generator.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace opensn;

namespace
{

constexpr unsigned int NUM_GROUPS = 2;

// The extruded mesh of transport_3d_4_cycles_1, whose cell graph has local cycles and, on the
// 2 x 2 KBA partition of 4 ranks, cycles between ranks.
std::shared_ptr<Mesh>
MakeCyclicMesh()
{
  ParameterBlock file_block;
  file_block.AddParameter(
    "filename", std::string(OPENSN_TEST_ROOT) + "/assets/mesh/square2x2_partition_cyclic3.obj");
  auto file_params = FromFileMeshGenerator::GetInputParameters();
  file_params.AssignParameters(file_block);
  const std::shared_ptr<MeshGenerator> from_file =
    std::make_shared<FromFileMeshGenerator>(file_params);

  ParameterBlock inputs("inputs");
  inputs.ChangeToArray();
  inputs.AddParameter(ParameterBlock("0", from_file));

  ParameterBlock layers("layers");
  layers.ChangeToArray();
  for (int i = 0; i < 4; ++i)
  {
    ParameterBlock layer(std::to_string(i));
    layer.AddParameter("z", 0.4 * (i + 1));
    layer.AddParameter("n", 2);
    layers.AddParameter(layer);
  }

  ParameterBlock block;
  block.AddParameter(inputs);
  block.AddParameter(layers);
  if (mpi_comm.size() == 4)
  {
    ParameterBlock kba_block;
    kba_block.AddParameter("nx", 2);
    kba_block.AddParameter("ny", 2);
    kba_block.AddParameter("xcuts", std::vector<double>{0.0});
    kba_block.AddParameter("ycuts", std::vector<double>{0.0});
    auto kba_params = KBAGraphPartitioner::GetInputParameters();
    kba_params.AssignParameters(kba_block);
    block.AddParameter(
      "partitioner",
      std::shared_ptr<GraphPartitioner>(std::make_shared<KBAGraphPartitioner>(kba_params)));
  }

  auto params = ExtruderMeshGenerator::GetInputParameters();
  params.AssignParameters(block);
  auto grid = ExtruderMeshGenerator(params).Execute();
  grid->SetUniformBlockID(0);
  return grid;
}

// Two-group absorber with in-group scattering.
std::shared_ptr<MultiGroupXS>
MakeXS()
{
  return std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromOpenSn(
    std::string(OPENSN_TEST_ROOT) + "/assets/xs/cbc_correction_2g.xs"));
}

// CBC problem with one two-group groupset and polar angle aggregation, so every angle set holds
// two directions.
std::shared_ptr<DiscreteOrdinatesProblem>
MakeProblem()
{
  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, NUM_GROUPS - 1}));
  groupset.AddParameter(
    "angular_quadrature",
    std::shared_ptr<AngularQuadrature>(std::make_shared<GLCProductQuadrature3DXYZ>(4, 8, 0)));
  groupset.AddParameter("angle_aggregation_type", std::string("polar"));
  groupset.AddParameter("inner_linear_method", std::string("petsc_gmres"));
  ParameterBlock groupsets("groupsets");
  groupsets.ChangeToArray();
  groupsets.AddParameter(groupset);

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{0}));
  xs_entry.AddParameter("xs", MakeXS());
  ParameterBlock xs_map("xs_map");
  xs_map.ChangeToArray();
  xs_map.AddParameter(xs_entry);

  ParameterBlock params;
  params.AddParameter("mesh", MakeCyclicMesh());
  params.AddParameter("num_groups", static_cast<int>(NUM_GROUPS));
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter("sweep_type", std::string("CBC"));
  return DiscreteOrdinatesProblem::Create(params);
}

// Distinct positive value for every (cell, face, downwind face node, group, direction).
double
Encode(std::uint32_t cell_local_id,
       unsigned int face,
       unsigned int face_node,
       unsigned int group,
       std::size_t angle)
{
  return 1.0 + group +
         NUM_GROUPS * (static_cast<double>(angle) +
                       256.0 * (face_node + 16.0 * (face + 16.0 * cell_local_id)));
}

} // namespace

// Verify that each nodal correction reaches the kernel's upwind entry for every
// delayed face, direction and group.
TEST(CBCFLUDSLaggedFluxCorrection, CorrectionLandsWhereTheSweepReadsIt)
{
  const auto problem = MakeProblem();
  const auto& grid = *problem->GetGrid();
  auto& angle_agg = *problem->GetGroupsets().front().angle_agg;

  const auto num_delayed = angle_agg.GetNewDelayedAngularDOFsAsSTLVector().size();
  std::size_t num_read = 0;
  std::size_t num_local_read = 0;
  std::size_t num_nonlocal_read = 0;
  std::size_t num_permuted_nodes = 0;
  std::size_t num_wrong = 0;
  double sum_read = 0.0;
  const auto VisitDelayed = [&](const auto& visit)
  {
    for (auto& angle_set : angle_agg)
    {
      auto& fluds = dynamic_cast<CBC_FLUDS&>(angle_set->GetFLUDS());
      const auto& common_data = fluds.GetCommonData();
      const auto num_angles = angle_set->GetAngleIndices().size();
      ASSERT_EQ(angle_set->GetNumGroups(), NUM_GROUPS);
      ASSERT_GT(num_angles, 1U);
      for (const auto& cell_ptr : grid.GetLocalCells())
      {
        const auto cell_local_id = static_cast<std::uint32_t>(cell_ptr->local_id);
        for (unsigned int f = 0; f < cell_ptr->faces.size(); ++f)
        {
          const bool local = common_data.IsDelayedLocalIncomingFace(cell_local_id, f);
          const bool nonlocal = common_data.IsDelayedNonlocalIncomingFace(cell_local_id, f);
          if (not local and not nonlocal)
            continue;
          const auto& mapping =
            common_data.GetFaceNodalMapping(cell_local_id, f).face_node_mapping_;
          for (unsigned int fj = 0; fj < mapping.size(); ++fj)
          {
            const auto upwind_node = static_cast<unsigned int>(mapping[fj]);
            if (upwind_node != fj)
              ++num_permuted_nodes;
            for (std::size_t n = 0; n < num_angles; ++n)
            {
              double* psi = local
                              ? fluds.DelayedUpwindPsi(cell_local_id, f, upwind_node, n)
                              : fluds.DelayedNLUpwindPsi(
                                  common_data.GetDelayedNonlocalFaceByLocalFace(cell_local_id, f),
                                  upwind_node,
                                  n);
              for (unsigned int g = 0; g < NUM_GROUPS; ++g)
              {
                visit(&psi[g], cell_local_id, f, fj, g, angle_set->GetAngleIndices()[n]);
                (local ? num_local_read : num_nonlocal_read) += 1;
              }
            }
          }
        }
      }
    }
  };

  std::unordered_map<const double*, double> initial;
  VisitDelayed(
    [&](double* psi, auto cell, auto face, auto node, auto group, auto angle)
    {
      const double value = Encode(cell, face, node, group, angle);
      EXPECT_TRUE(initial.emplace(psi, value).second);
      *psi = value;
    });
  EXPECT_EQ(initial.size(), num_delayed);
  angle_agg.SetDelayedPsiOld2New();
  // The correction API supplies one nodal delta to every direction.
  angle_agg.CorrectNewDelayedAngularFlux([](auto cell, auto face, auto node, auto group)
                                         { return Encode(cell, face, node, group, 0); });
  const auto corrected = angle_agg.GetNewDelayedAngularDOFsAsSTLVector();
  angle_agg.SetDelayedPsiNew2Old();

  std::unordered_set<const double*> visited;
  VisitDelayed(
    [&](double* psi, auto cell, auto face, auto node, auto group, auto)
    {
      EXPECT_TRUE(visited.insert(psi).second);
      const double expected = initial.at(psi) + Encode(cell, face, node, group, 0);
      if (*psi != expected)
        ++num_wrong;
      sum_read += *psi;
      ++num_read;
    });
  EXPECT_EQ(visited.size(), num_delayed);
  EXPECT_EQ(num_wrong, 0U);

  // Each delayed unknown has a unique address, is read once, and preserves its initial value
  // plus the nodal correction. The values also add up to those read.
  EXPECT_EQ(num_read, corrected.size());
  double sum_corrected = 0.0;
  std::size_t num_uncorrected = 0;
  for (const double value : corrected)
  {
    sum_corrected += value;
    if (value <= 0.0)
      ++num_uncorrected;
  }
  EXPECT_EQ(num_uncorrected, 0U);
  EXPECT_EQ(sum_corrected, sum_read);

  // The problem must exercise lagged local faces, a non-identity face-node mapping and, on four
  // ranks, lagged faces between ranks.
  std::size_t global_local_read = 0;
  std::size_t global_nonlocal_read = 0;
  std::size_t global_permuted = 0;
  mpi_comm.all_reduce(num_local_read, global_local_read, mpi::op::sum<std::size_t>());
  mpi_comm.all_reduce(num_nonlocal_read, global_nonlocal_read, mpi::op::sum<std::size_t>());
  mpi_comm.all_reduce(num_permuted_nodes, global_permuted, mpi::op::sum<std::size_t>());
  EXPECT_GT(global_local_read, 0U);
  EXPECT_GT(global_permuted, 0U);
  if (mpi_comm.size() == 4)
    EXPECT_GT(global_nonlocal_read, 0U);
}
