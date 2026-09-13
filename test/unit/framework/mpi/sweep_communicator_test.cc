// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mpi/sweep_communicator.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <limits>

using namespace opensn;

TEST(SweepCommunicatorTest, RankAndTagValidation)
{
  SweepCommunicator communicator(mpi_comm);

  EXPECT_EQ(communicator.GetPeerRank(mpi_comm.rank()), mpi_comm.rank());
  EXPECT_THROW(communicator.GetPeerRank(-1), std::out_of_range);
  EXPECT_THROW(communicator.GetPeerRank(mpi_comm.size()), std::out_of_range);

  EXPECT_EQ(communicator.BuildMessageTag(3, 8, 7), 31);
  EXPECT_THROW(communicator.BuildMessageTag(0, 0, 0), std::invalid_argument);
  EXPECT_THROW(communicator.BuildMessageTag(0, 2, 2), std::invalid_argument);
  EXPECT_THROW(communicator.BuildMessageTag(std::numeric_limits<std::size_t>::max(), 2, 1),
               std::out_of_range);
}

TEST(SweepCommunicatorTest, HasPrivateMessageContext)
{
  if (mpi_comm.size() < 2)
    return;

  SweepCommunicator sweep_communicator(mpi_comm);
  constexpr int tag = 7;

  mpi_comm.barrier();
  if (mpi_comm.rank() == 0)
  {
    const int world_value = 11;
    const int sweep_value = 22;
    auto world_request = mpi_comm.isend(1, tag, &world_value, 1);
    auto sweep_request = sweep_communicator.GetCommunicator().isend(1, tag, &sweep_value, 1);
    mpi::wait(world_request);
    mpi::wait(sweep_request);
  }
  else if (mpi_comm.rank() == 1)
  {
    int sweep_value = 0;
    int world_value = 0;
    sweep_communicator.GetCommunicator().recv(0, tag, &sweep_value, 1);
    mpi_comm.recv(0, tag, &world_value, 1);
    EXPECT_EQ(sweep_value, 22);
    EXPECT_EQ(world_value, 11);
  }
  mpi_comm.barrier();
}
