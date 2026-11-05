// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#ifdef OPENSN_WITH_PRECICE

#include "adapters/precice/cell_power.h"
#include "adapters/precice/precice_adapter.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/parameters/parameter_block.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/volumetric_source/volumetric_source.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"
#include <precice/precice.hpp>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

// These tests run OpenSn and a reference receiver as two preCICE participants in one MPI job. The
// first ranks form the OpenSn participant and the others the receiver; each side uses its own
// communicator as opensn::mpi_comm. They need at least two ranks and are skipped on one.

using namespace opensn;

namespace
{
using PointPowers = std::vector<double>;

// Reference receiver points along x for the synthetic sender.
const std::vector<Vector3> line_points{
  {0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};

/**
 * Receives power on `receiver_points`, distributed round-robin over the receiver ranks; a 2D
 * coupling mesh uses their x and y. With `writes_temperature`, also writes a constant temperature,
 * as an implicit scheme requires.
 */
PointPowers
RunReceiver(const std::filesystem::path& configuration,
            const std::vector<Vector3>& receiver_points,
            const bool writes_temperature)
{
  constexpr const char* mesh_name = "OpenFOAM-Cells";
  MPI_Comm communicator = mpi_comm;
  precice::Participant participant("OpenFOAM",
                                   configuration.string(),
                                   mpi_comm.rank(),
                                   mpi_comm.size(),
                                   static_cast<void*>(&communicator));
  const int dimension = participant.getMeshDimensions(mesh_name);

  std::vector<std::size_t> points;
  std::vector<double> coordinates;
  for (std::size_t point = 0; point < receiver_points.size(); ++point)
    if (point % static_cast<std::size_t>(mpi_comm.size()) ==
        static_cast<std::size_t>(mpi_comm.rank()))
    {
      points.push_back(point);
      for (int d = 0; d < dimension; ++d)
        coordinates.push_back(receiver_points[point][d]);
    }
  std::vector<precice::VertexID> vertex_ids(points.size());
  participant.setMeshVertices(mesh_name, coordinates, vertex_ids);
  participant.initialize();

  std::vector<double> powers(points.size(), 0.0);
  const std::vector<double> temperatures(points.size(), 300.0);
  while (participant.isCouplingOngoing())
  {
    participant.requiresWritingCheckpoint();
    const double time_step = participant.getMaxTimeStepSize();
    participant.readData(mesh_name, "TransportCellPower", vertex_ids, time_step, powers);
    if (writes_temperature)
      participant.writeData(mesh_name, "Temperature", vertex_ids, temperatures);
    participant.advance(time_step);
    participant.requiresReadingCheckpoint();
  }
  participant.finalize();

  PointPowers received(receiver_points.size(), 0.0);
  for (std::size_t i = 0; i < points.size(); ++i)
    received[points[i]] = powers[i];
  return received;
}

/**
 * A private directory holding a copy of a test configuration, which is every rank's working
 * directory while this object exists. preCICE writes its connection, iteration, and profiling
 * files to the working directory, so stale or concurrent runs cannot interfere and nothing is left
 * behind. Construction and destruction are collective.
 */
class RunDirectory
{
public:
  explicit RunDirectory(const std::string& configuration_name)
    : previous_(std::filesystem::current_path())
  {
    std::string directory;
    if (mpi_comm.rank() == 0)
    {
      const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
      directory =
        (std::filesystem::temp_directory_path() / ("opensn_precice_" + std::to_string(stamp)))
          .string();
      std::filesystem::create_directories(directory);

      std::ifstream input(std::filesystem::path(OPENSN_TEST_ROOT) / "unit" / "adapters" /
                          "precice" / configuration_name);
      std::stringstream buffer;
      buffer << input.rdbuf();
      std::string configuration = buffer.str();
      // Log only warnings and errors.
      const std::string root = "<precice-configuration>";
      const auto position = configuration.find(root);
      if (position == std::string::npos)
        throw std::runtime_error("test configuration has no precice-configuration element");
      configuration.insert(position + root.size(),
                           "\n  <log>\n    <sink filter=\"%Severity% >= warning\" "
                           "output=\"stdout\" />\n  </log>");
      std::ofstream(std::filesystem::path(directory) / "precice-config.xml") << configuration;
    }
    int length = static_cast<int>(directory.size());
    mpi_comm.broadcast(length, 0);
    directory.resize(length);
    mpi_comm.broadcast(directory.data(), length, 0);
    directory_ = directory;
    std::filesystem::current_path(directory_);
  }

  ~RunDirectory()
  {
    std::filesystem::current_path(previous_);
    mpi_comm.barrier();
    if (mpi_comm.rank() == 0)
      std::filesystem::remove_all(directory_);
  }

  RunDirectory(const RunDirectory&) = delete;
  RunDirectory& operator=(const RunDirectory&) = delete;

  std::filesystem::path GetConfiguration() const { return directory_ / "precice-config.xml"; }

private:
  std::filesystem::path previous_;
  std::filesystem::path directory_;
};

/**
 * Runs `sender` as the OpenSn participant on the first `num_senders` ranks and the reference
 * receiver on the others, with the named configuration from this directory. Returns the power
 * received at each receiver point on every rank.
 */
PointPowers
RunCoupled(const std::string& configuration_name,
           const int num_senders,
           const std::function<void(const std::filesystem::path&)>& sender,
           const std::vector<Vector3>& receiver_points,
           const bool receiver_writes_temperature = false)
{
  const mpi::Communicator world = mpi_comm;
  const int rank = world.rank();
  const RunDirectory run_directory(configuration_name);
  const auto configuration = run_directory.GetConfiguration();

  const bool is_sender = rank < num_senders;
  MPI_Comm participant_communicator = MPI_COMM_NULL;
  MPI_Comm_split(world, is_sender ? 0 : 1, rank, &participant_communicator);
  PointPowers local(receiver_points.size(), 0.0);
  mpi_comm = participant_communicator;
  try
  {
    if (is_sender)
      sender(configuration);
    else
      local = RunReceiver(configuration, receiver_points, receiver_writes_temperature);
  }
  catch (const std::exception& error)
  {
    // The other participant would wait for this one indefinitely.
    std::cerr << "Rank " << rank << ": " << error.what() << std::endl;
    world.abort(EXIT_FAILURE);
  }
  mpi_comm = world;
  MPI_Comm_free(&participant_communicator);

  PointPowers received(receiver_points.size(), 0.0);
  world.all_reduce(local.data(),
                   static_cast<int>(receiver_points.size()),
                   received.data(),
                   mpi::op::sum<double>());
  return received;
}

int
NumSenders()
{
  return (mpi_comm.size() + 1) / 2;
}

/// Powers 1, 2, and 3 W at x = 0, 1, 2 on all but the last sender rank, which owns no cells.
coupling::CellPowerData
MakeSyntheticPower(const unsigned int dimension)
{
  const auto owners = static_cast<std::uint64_t>(std::max(1, mpi_comm.size() - 1));
  std::vector<coupling::OwnedCellPower> cells;
  for (std::uint64_t cell = 0; cell < 3; ++cell)
    if (cell * owners / 3 == static_cast<std::uint64_t>(mpi_comm.rank()))
      cells.push_back({cell,
                       0,
                       Vector3(static_cast<double>(cell), 0.0, 0.0),
                       1.0,
                       static_cast<double>(cell + 1)});
  return {dimension, std::move(cells)};
}

void
ExpectPowers(const PointPowers& received, const PointPowers& expected)
{
  ASSERT_EQ(received.size(), expected.size());
  for (std::size_t point = 0; point < expected.size(); ++point)
    EXPECT_NEAR(received[point], expected[point], 1.0e-12) << "receiver point " << point;
}

// Data of the synthetic three-group library test/assets/xs/cepxs_synthetic_csda_3g.bxslib, from
// test/assets/xs/generate_cepxs_synthetic.py: total cross section (cm^-1), energy-deposition cross
// section (MeV/cm), and the downscatter transfer matrix S0[to][from] (cm^-1).
constexpr double synthetic_sigma_t[3] = {1.0, 1.2, 1.5};
constexpr double synthetic_sigma_edep[3] = {0.5, 0.35, 0.2};
constexpr double synthetic_s0[3][3] = {{0.20, 0.00, 0.00}, {0.10, 0.15, 0.00}, {0.03, 0.12, 0.10}};

/**
 * Exact energy deposition (MeV/cm^3) of the synthetic library in an infinite medium with a unit
 * group-0 source and CSDA disabled. The flux solves (diag(sigma_t) - S0) phi = Q, which is lower
 * triangular, and the deposition is sum_g sigma_edep[g] phi[g].
 */
double
ExactInfiniteMediumDeposition()
{
  const double source[3] = {1.0, 0.0, 0.0};
  double phi[3] = {};
  double deposition = 0.0;
  for (int g = 0; g < 3; ++g)
  {
    double inscatter = source[g];
    for (int gp = 0; gp < g; ++gp)
      inscatter += synthetic_s0[g][gp] * phi[gp];
    phi[g] = inscatter / (synthetic_sigma_t[g] - synthetic_s0[g][g]);
    deposition += synthetic_sigma_edep[g] * phi[g];
  }
  return deposition;
}

/**
 * Solves the synthetic three-group library, with CSDA disabled, on [0, 4]^dim with reflecting
 * boundaries and a unit group-0 source, and returns the energy-deposition field function, which is
 * exactly ExactInfiniteMediumDeposition() everywhere. Block 1 is the half with first coordinate
 * below 2 (z in 1D, x otherwise); block 2 is the rest.
 */
std::shared_ptr<FieldFunctionGridBased>
SolveInfiniteMediumDeposition(const unsigned int dimension)
{
  auto grid = dimension == 1   ? BuildLineMesh(4.0, 4, 0.0)
              : dimension == 2 ? BuildSquareMesh(4.0, 4, 0.0)
                               : BuildBoxMesh(4.0, 4, 0.0);
  const auto block_of = [dimension](const Vector3& centroid)
  { return (dimension == 1 ? centroid.z : centroid.x) < 2.0 ? 1U : 2U; };
  for (const auto& cell : grid->GetLocalCells())
    cell->block_id = block_of(cell->centroid);
  for (const auto ghost_id : grid->GetGhostGlobalIDs())
  {
    auto& cell = grid->GetGlobalCell(ghost_id);
    cell.block_id = block_of(cell.centroid);
  }

  std::shared_ptr<AngularQuadrature> quadrature;
  std::vector<std::string> boundaries;
  if (dimension == 1)
  {
    quadrature = std::make_shared<GLProductQuadrature1DSlab>(4, 0);
    boundaries = {"zmin", "zmax"};
  }
  else if (dimension == 2)
  {
    quadrature = std::make_shared<GLCProductQuadrature2DXY>(2, 4, 0);
    boundaries = {"xmin", "xmax", "ymin", "ymax"};
  }
  else
  {
    quadrature = std::make_shared<GLCProductQuadrature3DXYZ>(2, 4, 0);
    boundaries = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
  }

  ParameterBlock groupset("0");
  groupset.AddParameter(ParameterBlock("groups_from_to", std::vector<int>{0, 2}));
  groupset.AddParameter("angular_quadrature", quadrature);
  groupset.AddParameter("inner_linear_method", std::string("petsc_gmres"));
  groupset.AddParameter("l_abs_tol", 1.0e-12);
  groupset.AddParameter("l_max_its", 200);
  ParameterBlock groupsets("groupsets");
  groupsets.ChangeToArray();
  groupsets.AddParameter(groupset);

  ParameterBlock xs_entry("0");
  xs_entry.AddParameter(ParameterBlock("block_ids", std::vector<int>{1, 2}));
  const auto library =
    std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "xs" / "cepxs_synthetic_csda_3g.bxslib";
  xs_entry.AddParameter(
    "xs", std::make_shared<MultiGroupXS>(MultiGroupXS::LoadFromCEPXS(library.string(), 0, true)));
  ParameterBlock xs_map("xs_map");
  xs_map.ChangeToArray();
  xs_map.AddParameter(xs_entry);

  ParameterBlock source_params;
  source_params.AddParameter(ParameterBlock("block_ids", std::vector<int>{1, 2}));
  source_params.AddParameter(ParameterBlock("group_strength", std::vector<double>{1.0, 0.0, 0.0}));
  ParameterBlock sources("volumetric_sources");
  sources.ChangeToArray();
  sources.AddParameter("0", VolumetricSource::Create(source_params));

  ParameterBlock boundary_conditions("boundary_conditions");
  boundary_conditions.ChangeToArray();
  for (std::size_t b = 0; b < boundaries.size(); ++b)
  {
    ParameterBlock boundary(std::to_string(b));
    boundary.AddParameter("name", boundaries[b]);
    boundary.AddParameter("type", std::string("reflecting"));
    boundary_conditions.AddParameter(boundary);
  }

  ParameterBlock params;
  params.AddParameter("mesh", grid);
  params.AddParameter("num_groups", 3);
  params.AddParameter(groupsets);
  params.AddParameter(xs_map);
  params.AddParameter(sources);
  params.AddParameter(boundary_conditions);
  auto problem = DiscreteOrdinatesProblem::Create(params);

  ParameterBlock solver_params;
  solver_params.AddParameter("problem", std::shared_ptr<Problem>(problem));
  auto solver = SteadyStateSourceSolver::Create(solver_params);
  solver->Initialize();
  solver->Execute();
  return problem->CreateFieldFunction("deposition", "energy_deposition", -1.0);
}
} // namespace

// Conservative nearest-neighbor mapping sends each cell's power to the nearest receiver point:
// x = 0, 1, 2 land on points 0, 2, and 3. With more than two ranks, one OpenSn rank owns no cells
// and the receiver is parallel. The adapter may be used once.
TEST(PreciceExchange, ConservesPowerAcrossParticipants)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "two participants need at least two ranks";
  const auto received = RunCoupled(
    "precice-config.xml",
    NumSenders(),
    [](const std::filesystem::path& configuration)
    {
      const auto cell_data = MakeSyntheticPower(3);
      coupling::PreciceAdapter adapter(configuration);
      adapter.WriteSteadyCellPower(cell_data);
      EXPECT_THROW(adapter.WriteSteadyCellPower(cell_data), std::logic_error);
    },
    line_points);
  ExpectPowers(received, {1.0, 0.0, 2.0, 3.0});
}

// Implicit schemes require the checkpoint queries every iteration and data in both directions;
// OpenSn ignores the temperature it receives.
TEST(PreciceExchange, SupportsImplicitCoupling)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "two participants need at least two ranks";
  const auto received = RunCoupled(
    "precice-config-implicit.xml",
    NumSenders(),
    [](const std::filesystem::path& configuration)
    {
      coupling::PreciceAdapter adapter(configuration);
      adapter.WriteSteadyCellPower(MakeSyntheticPower(3));
    },
    line_points,
    true);
  ExpectPowers(received, {1.0, 0.0, 2.0, 3.0});
}

// A 2D OpenSn mesh may use a 2D coupling mesh. A length scale of 0.5 moves the centroids to
// x = 0, 0.5, 1, whose nearest receiver points are 0, 1, and 2; unscaled, they would be 0, 2, 3.
TEST(PreciceExchange, ScalesCoordinatesOn2DCouplingMesh)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "two participants need at least two ranks";
  const auto received = RunCoupled(
    "precice-config-2d.xml",
    NumSenders(),
    [](const std::filesystem::path& configuration)
    {
      coupling::PreciceAdapter adapter(
        configuration, "OpenSn", "OpenSn-Cells", "TransportCellPower", 0.5);
      adapter.WriteSteadyCellPower(MakeSyntheticPower(2));
    },
    line_points);
  ExpectPowers(received, {1.0, 2.0, 3.0, 0.0});
}

// A 3D OpenSn mesh cannot use a 2D coupling mesh. The check fails before any connection, so the
// receiver does not participate.
TEST(PreciceExchange, Rejects3DCellsOn2DCouplingMesh)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "two participants need at least two ranks";
  const mpi::Communicator world = mpi_comm;
  const RunDirectory run_directory("precice-config-2d.xml");
  MPI_Comm single = MPI_COMM_NULL;
  MPI_Comm_split(world, world.rank() == 0 ? 0 : MPI_UNDEFINED, 0, &single);
  if (world.rank() == 0)
  {
    mpi_comm = single;
    {
      coupling::PreciceAdapter adapter(run_directory.GetConfiguration());
      EXPECT_THROW(adapter.WriteSteadyCellPower(MakeSyntheticPower(3)), std::invalid_argument);
    }
    mpi_comm = world;
    MPI_Comm_free(&single);
  }
}

// Transport through the full path: the energy-deposition field of a multigroup solve, block
// selection, the transverse measure, conservative group averaging, the length scale, and the
// exchange. Coupled power is source_rate * joules_per_energy_unit * q * V = 1.5 * q * V, where q is
// the exact infinite-medium deposition and V is the physical volume of block 1. The 1e-8 tolerance
// covers the 1e-12 transport iteration tolerance. The receiver points
// are cell centers of a 3-per-direction grid over block 1 after scaling by 0.5, so each OpenSn
// rank's coupled cells overlap them, as preCICE's partition filtering requires.
namespace
{
class PreciceTransportExchange : public ::testing::TestWithParam<unsigned int>
{
};
} // namespace

TEST_P(PreciceTransportExchange, DeliversExactInfiniteMediumDeposition)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "two participants need at least two ranks";
  const unsigned int dimension = GetParam();
  // 1D: 2 x 3 cross-sectional area; 2D: 8 x 0.5 depth; 3D: 32.
  const double transverse = dimension == 1 ? 3.0 : 0.5;
  const double coupled_volume = dimension == 1 ? 6.0 : dimension == 2 ? 4.0 : 32.0;
  // Scaled block 1 spans [0, 1] in its first coordinate and [0, 2] in the others.
  const std::vector<double> first{1.0 / 6.0, 0.5, 5.0 / 6.0};
  const std::vector<double> other{1.0 / 3.0, 1.0, 5.0 / 3.0};
  std::vector<Vector3> receiver_points;
  if (dimension == 1)
    for (const double z : first)
      receiver_points.emplace_back(0.0, 0.0, z);
  else
    for (const double x : first)
      for (const double y : other)
        for (const double z : dimension == 2 ? std::vector<double>{0.0} : other)
          receiver_points.emplace_back(x, y, z);
  const auto received = RunCoupled(
    dimension == 2 ? "precice-config-2d.xml" : "precice-config.xml",
    NumSenders(),
    [dimension, transverse](const std::filesystem::path& configuration)
    {
      const auto field = SolveInfiniteMediumDeposition(dimension);
      auto cell_data = coupling::MakeCellPowerData(
        *field, 0, 3.0, 0.5, dimension < 3 ? std::optional<double>(transverse) : std::nullopt, {1});
      std::vector<std::int64_t> groups;
      for (const auto& cell : cell_data.GetCells())
        groups.push_back((dimension == 1 ? cell.centroid.z : cell.centroid.x) < 1.0 ? 0 : 1);
      cell_data.ConservativelyAverageByGroup(groups);
      coupling::PreciceAdapter adapter(
        configuration, "OpenSn", "OpenSn-Cells", "TransportCellPower", 0.5);
      adapter.WriteSteadyCellPower(cell_data);
    },
    receiver_points);

  double total = 0.0;
  for (const double power : received)
    total += power;
  const double expected = 1.5 * ExactInfiniteMediumDeposition() * coupled_volume;
  EXPECT_NEAR(total, expected, 1.0e-8 * expected);
}

INSTANTIATE_TEST_SUITE_P(Dimensions, PreciceTransportExchange, ::testing::Values(1U, 2U, 3U));

#endif // OPENSN_WITH_PRECICE
