// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{
constexpr unsigned int FACE_VALUES = 8;
constexpr int WHOLE_MESSAGE_BYTES = 1024;
// One kind, slot, total, offset and count, with room for one double per fragment.
constexpr int FRAGMENT_MESSAGE_BYTES =
  sizeof(std::uint8_t) + 4 * sizeof(std::size_t) + sizeof(double);

// Test record rejection locally through the unpack operation, without aborting MPI.
class TestCommunicator : public CBC_AsynchronousCommunicator
{
public:
  using CBC_AsynchronousCommunicator::CBC_AsynchronousCommunicator;
  using CBC_AsynchronousCommunicator::UnpackBuffer;
};

template <typename T>
void
Put(std::vector<char>& message, T value)
{
  const auto size = message.size();
  message.resize(size + sizeof(T));
  std::memcpy(message.data() + size, &value, sizeof(T));
}

std::vector<char>
Record(CBC_MessageKind kind,
       std::size_t slot = 0,
       std::size_t total = FACE_VALUES,
       std::size_t offset = 0,
       std::size_t count = FACE_VALUES)
{
  std::vector<char> message;
  Put(message, static_cast<std::uint8_t>(kind));
  if (kind == CBC_MessageKind::DELAYED_COMPLETION)
    return message;
  Put(message, slot);
  if (kind == CBC_MessageKind::NORMAL_FACE_PSI_CHUNK or
      kind == CBC_MessageKind::DELAYED_FACE_PSI_CHUNK)
  {
    Put(message, total);
    Put(message, offset);
    Put(message, count);
  }
  for (std::size_t i = 0; i < count; ++i)
    Put(message, static_cast<double>(i + offset + 1));
  return message;
}

class CBCAsyncComm : public testing::Test
{
protected:
  void SetUp() override
  {
    if (mpi_comm.size() < 2)
      GTEST_SKIP() << "The face-stream fixture requires at least two MPI ranks";
    grid_ = Mesh::New();
    grid_->AddGlobalVertex(0, Vector3{0.0, 0.0, 0.0});
    grid_->AddGlobalVertex(1, Vector3{0.0, 0.0, 1.0});
    // Two reciprocal faces form a rank cycle. Cut 0 -> 1: rank 0 receives an ordinary face
    // from rank 1, rank 1 receives a delayed face from rank 0, and further ranks own no cells.
    for (unsigned int rank = 0; rank < 2; ++rank)
    {
      auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
      cell->global_id = rank;
      cell->partition_id = static_cast<int>(rank);
      cell->vertex_ids = {0, 1};
      for (unsigned int f = 0; f < 2; ++f)
      {
        CellFace face;
        face.normal = Vector3{0.0, 0.0, rank == f ? -1.0 : 1.0};
        face.vertex_ids = {f};
        face.has_neighbor = true;
        face.neighbor_id = 1 - rank;
        cell->faces.push_back(face);
      }
      grid_->AddGlobalCell(cell);
    }
    const auto face_info = BuildSPDSFaceNeighborInfo(*grid_);
    spds_ = std::make_unique<CBC_SPDS>(0, Vector3{0.0, 0.0, 1.0}, grid_, face_info, true);
    spds_->SetGlobalSweepFAS({0, 1});
    spds_->ApplyGlobalSweepFAS();
    // FaceNodalMapping has const members, so build the vector rather than assign its elements.
    mappings_ = std::vector<CellFaceNodalMapping>(
      grid_->GetLocalCellCount(),
      CellFaceNodalMapping{FaceNodalMapping(0, {0}, {0}), FaceNodalMapping(1, {0}, {1})});
    common_ = CBC_FLUDSCommonData::MakeAlpha(*spds_, mappings_);
    common_->FinalizeBeta();
    unknown_.AddUnknown(UnknownType::VECTOR_N, FACE_VALUES);
    sdm_ = PieceWiseLinearDiscontinuous::New(grid_);
    fluds_ = std::make_unique<CBC_FLUDS>(FACE_VALUES, 1, *common_, unknown_, *sdm_);
    sweep_comm_ = std::make_unique<SweepCommunicator>(mpi_comm, 1);
    comm_ = std::make_unique<TestCommunicator>(0, 0, *fluds_, WHOLE_MESSAGE_BYTES, *sweep_comm_);
    comm_->InitializeDelayedUpstreamData();
  }

  bool Active() const { return mpi_comm.rank() < 2; }
  bool Delayed() const { return mpi_comm.rank() == 1; }
  int Source() const { return 1 - mpi_comm.rank(); }
  CBC_MessageKind Whole() const
  {
    return Delayed() ? CBC_MessageKind::DELAYED_FACE_PSI : CBC_MessageKind::NORMAL_FACE_PSI;
  }
  CBC_MessageKind Chunk() const
  {
    return Delayed() ? CBC_MessageKind::DELAYED_FACE_PSI_CHUNK
                     : CBC_MessageKind::NORMAL_FACE_PSI_CHUNK;
  }
  void Reset()
  {
    comm_->Reset();
    fluds_->ClearLocalAndReceivePsi();
    cells_.clear();
  }
  void CheckValues()
  {
    if (not Active())
      return;
    if (Delayed())
      fluds_->SetDelayedOutgoingPsiNewToOld();
    const auto* psi = Delayed()
                        ? fluds_->DelayedNLUpwindPsi(common_->GetDelayedNonlocalFaceBySlot(0), 0, 0)
                        : fluds_->NLUpwindPsi(0, 0, 0);
    ASSERT_NE(psi, nullptr);
    for (unsigned int i = 0; i < FACE_VALUES; ++i)
      EXPECT_EQ(psi[i], static_cast<double>(i + 1));
  }

  std::shared_ptr<Mesh> grid_;
  std::vector<CellFaceNodalMapping> mappings_;
  std::unique_ptr<CBC_SPDS> spds_;
  std::unique_ptr<CBC_FLUDSCommonData> common_;
  UnknownManager unknown_;
  std::shared_ptr<PieceWiseLinearDiscontinuous> sdm_;
  std::unique_ptr<CBC_FLUDS> fluds_;
  std::unique_ptr<SweepCommunicator> sweep_comm_;
  std::unique_ptr<TestCommunicator> comm_;
  std::vector<std::uint32_t> cells_;
};
} // namespace

TEST_F(CBCAsyncComm, WholeFragmentsCompletionAndReset)
{
  if (not Active())
  {
    EXPECT_TRUE(comm_->ReceiveDelayedData());
    return;
  }
  for (int repeat = 0; repeat < 2; ++repeat)
  {
    Reset();
    EXPECT_NO_THROW(comm_->UnpackBuffer(Record(Whole()), Source(), Delayed(), cells_));
    CheckValues();
    EXPECT_EQ(cells_.size(), Delayed() ? 0U : 1U);
    Reset();
    EXPECT_NO_THROW(
      comm_->UnpackBuffer(Record(Chunk(), 0, FACE_VALUES, 0, 2), Source(), Delayed(), cells_));
    EXPECT_TRUE(cells_.empty());
    EXPECT_NO_THROW(comm_->UnpackBuffer(
      Record(Chunk(), 0, FACE_VALUES, 2, FACE_VALUES - 2), Source(), Delayed(), cells_));
    CheckValues();
    EXPECT_EQ(cells_.size(), Delayed() ? 0U : 1U);
    if (Delayed())
    {
      EXPECT_FALSE(comm_->ReceiveDelayedData());
      EXPECT_NO_THROW(
        comm_->UnpackBuffer(Record(CBC_MessageKind::DELAYED_COMPLETION), Source(), true, cells_));
    }
    EXPECT_TRUE(comm_->ReceiveDelayedData());
  }
}

TEST_F(CBCAsyncComm, EveryTruncatedHeaderAndPayload)
{
  if (not Active())
    return;
  for (const auto kind : {Whole(), Chunk()})
  {
    const auto message = Record(kind);
    for (std::size_t length = 0; length < message.size(); ++length)
    {
      Reset();
      EXPECT_THROW(
        comm_->UnpackBuffer(std::span(message).first(length), Source(), Delayed(), cells_),
        std::logic_error)
        << static_cast<int>(kind) << " length " << length;
      EXPECT_TRUE(cells_.empty());
    }
  }
}

TEST_F(CBCAsyncComm, InvalidKindsSourcesSlotsAndBounds)
{
  if (not Active())
    return;
  constexpr auto invalid_kind = static_cast<std::uint8_t>(CBC_MessageKind::DELAYED_COMPLETION) + 1;
  for (const auto& message :
       std::vector<std::vector<char>>{{invalid_kind},
                                      Record(Whole(), 1),
                                      Record(Whole(), std::numeric_limits<std::size_t>::max()),
                                      Record(Chunk(), 0, FACE_VALUES - 1),
                                      Record(Chunk(), 0, FACE_VALUES, 1, FACE_VALUES - 1),
                                      Record(Chunk(), 0, FACE_VALUES, 0, 0),
                                      Record(Chunk(), 0, std::numeric_limits<std::size_t>::max()),
                                      Record(Chunk(), 0, FACE_VALUES, 0, FACE_VALUES + 1)})
  {
    Reset();
    EXPECT_THROW(comm_->UnpackBuffer(message, Source(), Delayed(), cells_), std::logic_error);
    EXPECT_TRUE(cells_.empty());
  }
  EXPECT_THROW(comm_->UnpackBuffer(Record(Whole()), -1, Delayed(), cells_), std::logic_error);
  EXPECT_THROW(comm_->UnpackBuffer(Record(Whole()), mpi_comm.size(), Delayed(), cells_),
               std::logic_error);
  EXPECT_THROW(comm_->UnpackBuffer(Record(Whole()), mpi_comm.rank(), Delayed(), cells_),
               std::logic_error);
  EXPECT_THROW(comm_->UnpackBuffer(Record(Whole()), Source(), not Delayed(), cells_),
               std::logic_error);
  // Even when the kind matches the tag, a face cannot move to the other source/slot stream.
  const auto wrong_kind =
    Delayed() ? CBC_MessageKind::NORMAL_FACE_PSI : CBC_MessageKind::DELAYED_FACE_PSI;
  EXPECT_THROW(comm_->UnpackBuffer(Record(wrong_kind), Source(), not Delayed(), cells_),
               std::logic_error);
}

TEST_F(CBCAsyncComm, DuplicateOverlapAndCompletionMisuse)
{
  if (not Active())
    return;
  EXPECT_THROW(
    comm_->UnpackBuffer(Record(CBC_MessageKind::DELAYED_COMPLETION), Source(), Delayed(), cells_),
    std::logic_error);
  EXPECT_NO_THROW(
    comm_->UnpackBuffer(Record(Chunk(), 0, FACE_VALUES, 0, 2), Source(), Delayed(), cells_));
  EXPECT_THROW(
    comm_->UnpackBuffer(Record(Chunk(), 0, FACE_VALUES, 0, 2), Source(), Delayed(), cells_),
    std::logic_error);
  EXPECT_THROW(comm_->UnpackBuffer(
                 Record(Chunk(), 0, FACE_VALUES, 1, FACE_VALUES - 1), Source(), Delayed(), cells_),
               std::logic_error);
  EXPECT_NO_THROW(comm_->UnpackBuffer(
    Record(Chunk(), 0, FACE_VALUES, 2, FACE_VALUES - 2), Source(), Delayed(), cells_));
  EXPECT_THROW(comm_->UnpackBuffer(Record(Whole()), Source(), Delayed(), cells_), std::logic_error);
  if (Delayed())
  {
    EXPECT_NO_THROW(
      comm_->UnpackBuffer(Record(CBC_MessageKind::DELAYED_COMPLETION), Source(), true, cells_));
    EXPECT_THROW(
      comm_->UnpackBuffer(Record(CBC_MessageKind::DELAYED_COMPLETION), Source(), true, cells_),
      std::logic_error);
    EXPECT_THROW(comm_->UnpackBuffer(Record(Whole()), Source(), true, cells_), std::logic_error);
  }
}

TEST_F(CBCAsyncComm, PackedRecordsAndMalformedSuffix)
{
  if (not Active())
    return;
  auto message = Record(Whole());
  if (Delayed())
  {
    const auto completion = Record(CBC_MessageKind::DELAYED_COMPLETION);
    message.insert(message.end(), completion.begin(), completion.end());
    EXPECT_NO_THROW(comm_->UnpackBuffer(message, Source(), true, cells_));
    EXPECT_TRUE(comm_->ReceiveDelayedData());
  }
  else
  {
    // A later malformed record must be rejected without undoing the first stored face.
    message.push_back(static_cast<char>(Whole()));
    EXPECT_THROW(comm_->UnpackBuffer(message, Source(), false, cells_), std::logic_error);
    EXPECT_EQ(cells_, (std::vector<std::uint32_t>{0}));
  }
  CheckValues();
}

TEST_F(CBCAsyncComm, RejectZeroSizedFaceMetadata)
{
  // Zero directions make face payloads empty without dividing by zero groups in cell storage.
  CBC_FLUDS empty_values(FACE_VALUES, 0, *common_, unknown_, *sdm_);
  if (Active())
    EXPECT_THROW((TestCommunicator(0, 0, empty_values, WHOLE_MESSAGE_BYTES, *sweep_comm_)),
                 std::logic_error);
  else
    EXPECT_NO_THROW((TestCommunicator(0, 0, empty_values, WHOLE_MESSAGE_BYTES, *sweep_comm_)));
}

TEST_F(CBCAsyncComm, TransportWholeAndFragmentedRepeatedSweeps)
{
  // Extra ranks still construct, reset and drain an empty communicator.
  for (const int limit : {WHOLE_MESSAGE_BYTES, FRAGMENT_MESSAGE_BYTES})
  {
    comm_ = std::make_unique<TestCommunicator>(0, 0, *fluds_, limit, *sweep_comm_);
    comm_->InitializeDelayedUpstreamData();
    for (int repeat = 0; repeat < 2; ++repeat)
    {
      Reset();
      std::vector<double> values(FACE_VALUES);
      for (unsigned int i = 0; i < FACE_VALUES; ++i)
        values[i] = i + 1;
      if (Active())
      {
        const auto face = mpi_comm.rank() == 0 ? 1U : 0U;
        const bool delayed = mpi_comm.rank() == 0;
        comm_->QueueDownwindMessage(delayed ? TestCommunicator::DownwindPsiType::DELAYED
                                            : TestCommunicator::DownwindPsiType::NORMAL,
                                    common_->OutgoingPeerIndex(0, face),
                                    common_->OutgoingFaceSlot(0, face),
                                    values);
      }
      bool normal_received = mpi_comm.rank() != 0;
      std::size_t notifications = 0;
      for (;;)
      {
        comm_->ReceiveData(cells_);
        notifications += cells_.size();
        normal_received = normal_received or not cells_.empty();
        const bool sent = comm_->FlushSendBuffers();
        const bool delayed_received = comm_->ReceiveDelayedData();
        if (sent and delayed_received and normal_received)
          break;
      }
      EXPECT_EQ(notifications, mpi_comm.rank() == 0 ? 1U : 0U);
      CheckValues();
      // The next repetition reuses the same tags and must not consume another sweep's faces.
      mpi_comm.barrier();
    }
  }
}

TEST(CBCAsyncCommEmpty, ZeroWorkRank)
{
  const auto grid = Mesh::New();
  const auto face_info = BuildSPDSFaceNeighborInfo(*grid);
  CBC_SPDS spds(0, Vector3{0.0, 0.0, 1.0}, grid, face_info, true);
  const std::vector<CellFaceNodalMapping> mappings;
  auto common = CBC_FLUDSCommonData::MakeAlpha(spds, mappings);
  common->FinalizeBeta();
  const auto sdm = PieceWiseLinearDiscontinuous::New(grid);
  UnknownManager unknown;
  unknown.AddUnknown(UnknownType::VECTOR_N, FACE_VALUES);
  CBC_FLUDS fluds(FACE_VALUES, 1, *common, unknown, *sdm);
  SweepCommunicator sweep_comm(mpi_comm, 1);
  TestCommunicator comm(0, 0, fluds, WHOLE_MESSAGE_BYTES, sweep_comm);
  comm.InitializeDelayedUpstreamData();
  for (int repeat = 0; repeat < 2; ++repeat)
  {
    std::vector<std::uint32_t> cells;
    comm.ReceiveData(cells);
    EXPECT_TRUE(cells.empty());
    EXPECT_TRUE(comm.FlushSendBuffers());
    EXPECT_TRUE(comm.ReceiveDelayedData());
    EXPECT_THROW(
      comm.UnpackBuffer(Record(CBC_MessageKind::DELAYED_COMPLETION), mpi_comm.rank(), true, cells),
      std::logic_error);
    comm.Reset();
  }
}
