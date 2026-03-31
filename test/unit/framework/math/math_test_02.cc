// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gmock/gmock.h"
#include "framework/data_types/parallel_vector/ghosted_parallel_stl_vector.h"
#include "framework/data_types/parallel_vector/parallel_stl_vector.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

TEST(MathTest, ParallelVector)
{
  if (opensn::mpi_comm.size() != 2)
    return;

  ParallelSTLVector vec(5, 10, opensn::mpi_comm);

  if (opensn::mpi_comm.rank() == 0)
    vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    vec.SetValue(0, 1.0, VecOpType::SET_VALUE);
  vec.Assemble();

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_THAT(vec.GetLocalSTLData(), ::testing::ElementsAre(1, 0, 0, 0, 0));
  else
    EXPECT_THAT(vec.GetLocalSTLData(), ::testing::ElementsAre(2, 0, 0, 0, 0));

  const uint64_t ghost_id = opensn::mpi_comm.rank() == 0 ? 5 : 4;
  GhostedParallelSTLVector ghost_vec(5, 10, {ghost_id}, opensn::mpi_comm);

  EXPECT_EQ(ghost_vec.GetNumGhosts(), 1);

  if (opensn::mpi_comm.rank() == 0)
    ghost_vec.SetValue(5, 2.0, VecOpType::SET_VALUE);
  else
    ghost_vec.SetValue(4, 1.0, VecOpType::SET_VALUE);
  ghost_vec.Assemble();
  ghost_vec.CommunicateGhostEntries();

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_THAT(ghost_vec.GetLocalSTLData(), ::testing::ElementsAre(0, 0, 0, 0, 1, 2));
  else
    EXPECT_THAT(ghost_vec.GetLocalSTLData(), ::testing::ElementsAre(2, 0, 0, 0, 0, 1));

  {
    const auto made_vals = ghost_vec.MakeLocalVector();
    if (opensn::mpi_comm.rank() == 0)
      EXPECT_THAT(made_vals, ::testing::ElementsAre(0, 0, 0, 0, 1));
    else
      EXPECT_THAT(made_vals, ::testing::ElementsAre(2, 0, 0, 0, 0));
  }

  EXPECT_NEAR(vec.ComputeNorm(opensn::NormType::L1_NORM), 3., 1e-10);
  EXPECT_NEAR(vec.ComputeNorm(opensn::NormType::L2_NORM), 2.236067977, 1e-8);
  EXPECT_NEAR(vec.ComputeNorm(opensn::NormType::LINF_NORM), 2., 1e-10);

  EXPECT_NEAR(ghost_vec.ComputeNorm(opensn::NormType::L1_NORM), 3., 1e-10);
  EXPECT_NEAR(ghost_vec.ComputeNorm(opensn::NormType::L2_NORM), 2.236067977, 1e-8);
  EXPECT_NEAR(ghost_vec.ComputeNorm(opensn::NormType::LINF_NORM), 2., 1e-10);

  ParallelSTLVector vec2(5, 10, opensn::mpi_comm);

  vec2.CopyLocalValues(vec);

  if (opensn::mpi_comm.rank() == 0)
    vec2.SetValue(5, 2.0, VecOpType::ADD_VALUE);
  else
    vec2.SetValue(0, 1.0, VecOpType::ADD_VALUE);
  vec2.Assemble();

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_THAT(vec2.GetLocalSTLData(), ::testing::ElementsAre(2, 0, 0, 0, 0));
  else
    EXPECT_THAT(vec2.GetLocalSTLData(), ::testing::ElementsAre(4, 0, 0, 0, 0));

  ParallelSTLVector vec3(5, 10, opensn::mpi_comm);

  if (opensn::mpi_comm.rank() == 0)
    vec3.SetValues({5, 6}, {2.0, 3.0}, VecOpType::ADD_VALUE);
  else
    vec3.SetValues({0, 1}, {1.0, 4.0}, VecOpType::ADD_VALUE);
  vec3.Assemble();

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_THAT(vec3.GetLocalSTLData(), ::testing::ElementsAre(1, 4, 0, 0, 0));
  else
    EXPECT_THAT(vec3.GetLocalSTLData(), ::testing::ElementsAre(2, 3, 0, 0, 0));

  std::vector<uint64_t> ghost_ids;
  if (opensn::mpi_comm.rank() == 0)
    ghost_ids = {5, 6};
  else
    ghost_ids = {0, 1, 3};
  VectorGhostCommunicator vgc(5, 10, ghost_ids, opensn::mpi_comm);

  GhostedParallelSTLVector ghost_vec2(vgc);

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_EQ(ghost_vec2.GetLocalSizeWithGhosts(), 7);
  else
    EXPECT_EQ(ghost_vec2.GetLocalSizeWithGhosts(), 8);

  if (opensn::mpi_comm.rank() == 0)
    ghost_vec2.SetValues({5, 6}, {6.0, 7.0}, VecOpType::ADD_VALUE);
  else
    ghost_vec2.SetValues({0, 1, 3}, {1.0, 2.0, 4.0}, VecOpType::ADD_VALUE);

  ghost_vec2.Assemble();
  ghost_vec2.CommunicateGhostEntries();

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_THAT(ghost_vec2.GetLocalSTLData(), ::testing::ElementsAre(1, 2, 0, 4, 0, 6, 7));
  else
    EXPECT_THAT(ghost_vec2.GetLocalSTLData(), ::testing::ElementsAre(6, 7, 0, 0, 0, 1, 2, 4));

  {
    if (opensn::mpi_comm.rank() == 0)
      EXPECT_THAT(ghost_vec2.GetGhostIndices(), ::testing::ElementsAre(5, 6));
    else
      EXPECT_THAT(ghost_vec2.GetGhostIndices(), ::testing::ElementsAre(0, 1, 3));
  }

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_EQ(ghost_vec2.MapGhostToLocal(6), 6);
  else
    EXPECT_EQ(ghost_vec2.MapGhostToLocal(1), 6);

  {
    const auto ghosted_local = ghost_vec2.MakeGhostedLocalVector();
    if (opensn::mpi_comm.rank() == 0)
      EXPECT_THAT(ghosted_local, ::testing::ElementsAre(1, 2, 0, 4, 0, 6, 7));
    else
      EXPECT_THAT(ghosted_local, ::testing::ElementsAre(6, 7, 0, 0, 0, 1, 2, 4));
  }

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_NEAR(ghost_vec2.GetGlobalValue(3), 4., 1e-10);
  else
    EXPECT_NEAR(ghost_vec2.GetGlobalValue(6), 7., 1e-10);

  if (opensn::mpi_comm.rank() == 0)
    EXPECT_NEAR(ghost_vec2.GetGlobalValue(6), 7., 1e-10);
  else
    EXPECT_NEAR(ghost_vec2.GetGlobalValue(1), 2., 1e-10);
}
