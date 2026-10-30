#include "framework/graphs/kba_graph_partitioner.h"
#include <gmock/gmock.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <array>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

using namespace opensn;

TEST(KBAGraphPartitionerTest, Partition)
{
  ParameterBlock input_parameters;

  input_parameters.AddParameter("nx", 2);
  input_parameters.AddParameter("ny", 2);
  input_parameters.AddParameter("nz", 2);

  input_parameters.AddParameter("xcuts", std::vector<double>{0.0});
  input_parameters.AddParameter("ycuts", std::vector<double>{0.0});
  input_parameters.AddParameter("zcuts", std::vector<double>{0.0});

  InputParameters valid_parameters = KBAGraphPartitioner::GetInputParameters();

  valid_parameters.AssignParameters(input_parameters);

  KBAGraphPartitioner partitioner(valid_parameters);

  std::vector<std::vector<uint64_t>> dummy_graph(8);
  std::vector<Vector3> centroids = {{-1.0, -1.0, -1.0},
                                    {1.0, -1.0, -1.0},
                                    {-1.0, 1.0, -1.0},
                                    {1.0, 1.0, -1.0},
                                    {-1.0, -1.0, 1.0},
                                    {1.0, -1.0, 1.0},
                                    {-1.0, 1.0, 1.0},
                                    {1.0, 1.0, 1.0}};

  auto cell_pids = partitioner.Partition(dummy_graph, centroids, 2 * 2 * 2);
  EXPECT_THAT(cell_pids, testing::ElementsAre(0, 1, 2, 3, 4, 5, 6, 7));
}

namespace
{

KBAGraphPartitioner
MakeKBA(int nx, int ny, int nz)
{
  const auto cuts = [](int n)
  {
    std::vector<double> c;
    for (int i = 1; i < n; ++i)
      c.push_back(static_cast<double>(i) / n);
    return c;
  };
  ParameterBlock input_parameters;
  input_parameters.AddParameter("nx", nx);
  input_parameters.AddParameter("ny", ny);
  input_parameters.AddParameter("nz", nz);
  input_parameters.AddParameter("xcuts", cuts(nx));
  input_parameters.AddParameter("ycuts", cuts(ny));
  input_parameters.AddParameter("zcuts", cuts(nz));
  InputParameters valid_parameters = KBAGraphPartitioner::GetInputParameters();
  valid_parameters.AssignParameters(input_parameters);
  return KBAGraphPartitioner(valid_parameters);
}

/// Centroids of an n-by-n grid of unit-square cells in the z = 0 plane.
std::vector<Vector3>
PlanarCentroids(int n)
{
  std::vector<Vector3> centroids;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      centroids.emplace_back((i + 0.5) / n, (j + 0.5) / n, 0.0);
  return centroids;
}

} // namespace

TEST(KBAGraphPartitionerTest, PartitionsPlanarMeshInXAndY)
{
  auto partitioner = MakeKBA(2, 2, 1);
  const auto centroids = PlanarCentroids(2);
  const std::vector<std::vector<uint64_t>> graph(centroids.size());
  EXPECT_THAT(partitioner.Partition(graph, centroids, 4), testing::ElementsAre(0, 1, 2, 3));
}

TEST(KBAGraphPartitionerTest, RejectsZCutsOnPlanarMesh)
{
  auto partitioner = MakeKBA(2, 2, 2);
  const auto centroids = PlanarCentroids(4);
  const std::vector<std::vector<uint64_t>> graph(centroids.size());
  EXPECT_THROW(partitioner.Partition(graph, centroids, 8), std::invalid_argument);
}

TEST(KBAGraphPartitionerTest, RejectsXCutsOnLineMesh)
{
  // A 1D mesh lies along z.
  auto partitioner = MakeKBA(2, 1, 2);
  std::vector<Vector3> centroids;
  centroids.reserve(4);
  for (int k = 0; k < 4; ++k)
    centroids.emplace_back(0.0, 0.0, (k + 0.5) / 4.0);
  const std::vector<std::vector<uint64_t>> graph(centroids.size());
  EXPECT_THROW(partitioner.Partition(graph, centroids, 4), std::invalid_argument);
}

namespace
{

/// Partition of one centroid at the center of each cellset of an nx-by-ny-by-nz KBA cellset array
/// on the unit cube, listed in cellset order (x fastest).
std::vector<int>
PartitionCellsetCenters(int nx, int ny, int nz, int num_parts)
{
  auto partitioner = MakeKBA(nx, ny, nz);
  std::vector<Vector3> centroids;
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        centroids.emplace_back((i + 0.5) / nx, (j + 0.5) / ny, (k + 0.5) / nz);
  const std::vector<std::vector<uint64_t>> graph(centroids.size());
  return partitioner.Partition(graph, centroids, num_parts);
}

} // namespace

TEST(KBAGraphPartitionerTest, ConsecutiveCellsetGroupingIsKept)
{
  // Two cellsets per partition along x, as consecutive cellset ids give.
  EXPECT_THAT(PartitionCellsetCenters(2, 2, 2, 4), testing::ElementsAre(0, 0, 1, 1, 2, 2, 3, 3));
  // Whole rows per partition.
  EXPECT_THAT(PartitionCellsetCenters(2, 2, 1, 2), testing::ElementsAre(0, 0, 1, 1));
  // A serial run owns everything.
  EXPECT_THAT(PartitionCellsetCenters(2, 2, 2, 1), testing::ElementsAre(0, 0, 0, 0, 0, 0, 0, 0));
}

TEST(KBAGraphPartitionerTest, SubdomainsAreBricksOfCellsets)
{
  // With 3x2 cellsets on 3 partitions, consecutive ids would pair cellset (2, 0) with cellset
  // (0, 1). Each partition instead owns one column of two cellsets.
  EXPECT_THAT(PartitionCellsetCenters(3, 2, 1, 3), testing::ElementsAre(0, 1, 2, 0, 1, 2));
}

TEST(KBAGraphPartitionerTest, UngroupableCellsetsFallBackToConsecutiveIds)
{
  // Three cellsets cannot be split into two equal groups.
  EXPECT_THAT(PartitionCellsetCenters(3, 1, 1, 2), testing::ElementsAre(0, 0, 1));
}

namespace
{

KBAGraphPartitioner
MakeKBAWithCuts(const std::vector<double>& xcuts,
                const std::vector<double>& ycuts,
                const std::vector<double>& zcuts)
{
  ParameterBlock input_parameters;
  input_parameters.AddParameter("nx", static_cast<int>(xcuts.size()) + 1);
  input_parameters.AddParameter("ny", static_cast<int>(ycuts.size()) + 1);
  input_parameters.AddParameter("nz", static_cast<int>(zcuts.size()) + 1);
  input_parameters.AddParameter("xcuts", xcuts);
  input_parameters.AddParameter("ycuts", ycuts);
  input_parameters.AddParameter("zcuts", zcuts);
  InputParameters valid_parameters = KBAGraphPartitioner::GetInputParameters();
  valid_parameters.AssignParameters(input_parameters);
  return KBAGraphPartitioner(valid_parameters);
}

/// An nx-by-ny-by-nz grid of unit cells with their centroids and face adjacency. Cells are
/// numbered with x fastest.
struct CellGrid
{
  std::array<int, 3> n;
  std::vector<Vector3> centroids;
  std::vector<std::vector<uint64_t>> graph;

  uint64_t Index(int i, int j, int k) const
  {
    return static_cast<uint64_t>(i) +
           static_cast<uint64_t>(n[0]) *
             (static_cast<uint64_t>(j) + static_cast<uint64_t>(n[1]) * static_cast<uint64_t>(k));
  }
};

CellGrid
MakeCellGrid(int nx, int ny, int nz)
{
  CellGrid grid{{nx, ny, nz}, {}, {}};
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
      {
        grid.centroids.emplace_back(i + 0.5, j + 0.5, k + 0.5);
        std::vector<uint64_t> neighbors;
        const std::array<int, 3> ijk = {i, j, k};
        for (int d = 0; d < 3; ++d)
          for (const int step : {-1, 1})
          {
            auto nb = ijk;
            nb[d] += step;
            if (nb[d] >= 0 and nb[d] < grid.n[d])
              neighbors.push_back(grid.Index(nb[0], nb[1], nb[2]));
          }
        grid.graph.push_back(neighbors);
      }
  return grid;
}

/// True when the cells of every partition fill exactly their bounding brick of cells.
bool
PartitionsAreBricks(const CellGrid& grid, const std::vector<int>& pids, int num_parts)
{
  std::vector<std::array<int, 6>> bounds(num_parts,
                                         {std::numeric_limits<int>::max(),
                                          std::numeric_limits<int>::max(),
                                          std::numeric_limits<int>::max(),
                                          -1,
                                          -1,
                                          -1});
  std::vector<int> count(num_parts, 0);
  for (int k = 0; k < grid.n[2]; ++k)
    for (int j = 0; j < grid.n[1]; ++j)
      for (int i = 0; i < grid.n[0]; ++i)
      {
        const int p = pids[grid.Index(i, j, k)];
        const std::array<int, 3> ijk = {i, j, k};
        for (int d = 0; d < 3; ++d)
        {
          bounds[p][d] = std::min(bounds[p][d], ijk[d]);
          bounds[p][d + 3] = std::max(bounds[p][d + 3], ijk[d]);
        }
        ++count[p];
      }
  for (int p = 0; p < num_parts; ++p)
  {
    const auto& b = bounds[p];
    if (count[p] != (b[3] - b[0] + 1) * (b[4] - b[1] + 1) * (b[5] - b[2] + 1))
      return false;
  }
  return true;
}

/// True when, for every octant, the partition-level sweep dependency graph has no cycle.
bool
RankGraphAcyclicInEveryOctant(const CellGrid& grid, const std::vector<int>& pids, int num_parts)
{
  using RankGraph = boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>;
  for (int octant = 0; octant < 8; ++octant)
  {
    RankGraph g(num_parts);
    for (int k = 0; k < grid.n[2]; ++k)
      for (int j = 0; j < grid.n[1]; ++j)
        for (int i = 0; i < grid.n[0]; ++i)
        {
          const std::array<int, 3> ijk = {i, j, k};
          for (int d = 0; d < 3; ++d)
          {
            // Downstream neighbor along axis d for this octant's sign of Omega_d.
            auto nb = ijk;
            nb[d] += (octant & (1 << d)) != 0 ? 1 : -1;
            if (nb[d] < 0 or nb[d] >= grid.n[d])
              continue;
            const int up = pids[grid.Index(i, j, k)];
            const int down = pids[grid.Index(nb[0], nb[1], nb[2])];
            if (up != down)
              boost::add_edge(up, down, g);
          }
        }
    std::vector<RankGraph::vertex_descriptor> order;
    try
    {
      boost::topological_sort(g, std::back_inserter(order)); // NOLINT
    }
    catch (const boost::not_a_dag&)
    {
      return false;
    }
  }
  return true;
}

} // namespace

TEST(KBAGraphPartitionerTest, RejectsInvalidCuts)
{
  // A duplicate of the first cut, a decreasing cut, and a non-finite cut.
  EXPECT_THROW(MakeKBAWithCuts({0.5, 0.5}, {}, {}), std::invalid_argument);
  EXPECT_THROW(MakeKBAWithCuts({0.5, 0.25}, {}, {}), std::invalid_argument);
  EXPECT_THROW(MakeKBAWithCuts({std::numeric_limits<double>::quiet_NaN()}, {}, {}),
               std::invalid_argument);
  EXPECT_THROW(MakeKBA(0, 1, 1), std::invalid_argument);
}

TEST(KBAGraphPartitionerTest, RejectsSlabWithoutCells)
{
  // The second x cut lies beyond every cell centroid.
  const auto grid = MakeCellGrid(4, 2, 1);
  auto partitioner = MakeKBAWithCuts({2.0, 10.0}, {}, {});
  EXPECT_THROW(partitioner.Partition(grid.graph, grid.centroids, 3), std::invalid_argument);
}

TEST(KBAGraphPartitionerTest, PartitionsThinAxis)
{
  // A real extent is accepted however small it is next to the other axes: x spans 4e-10, y 500.
  const std::vector<Vector3> centroids = {
    {2.0e-10, 250.0, 0.0}, {6.0e-10, 250.0, 0.0}, {2.0e-10, 750.0, 0.0}, {6.0e-10, 750.0, 0.0}};
  const std::vector<std::vector<uint64_t>> graph(centroids.size());
  auto partitioner = MakeKBAWithCuts({4.0e-10}, {500.0}, {});
  EXPECT_THAT(partitioner.Partition(graph, centroids, 4), testing::ElementsAre(0, 1, 2, 3));
}

TEST(KBAGraphPartitionerTest, RejectsZCutsOnPlanarMeshOffZero)
{
  // Centroids of a triangle and a quadrilateral in the plane z = 0.1 average three and four equal
  // vertex coordinates, which round differently.
  const double z_triangle = (0.1 + 0.1 + 0.1) / 3.0;
  const double z_quadrilateral = (0.1 + 0.1 + 0.1 + 0.1) / 4.0;
  ASSERT_NE(z_triangle, z_quadrilateral);
  const std::vector<Vector3> centroids = {{0.25, 0.25, z_triangle},
                                          {0.75, 0.25, z_quadrilateral},
                                          {0.25, 0.75, z_quadrilateral},
                                          {0.75, 0.75, z_triangle}};
  const std::vector<std::vector<uint64_t>> graph(centroids.size());
  auto partitioner = MakeKBAWithCuts({0.5}, {0.5}, {0.1});
  EXPECT_THROW(partitioner.Partition(graph, centroids, 8), std::invalid_argument);
}

TEST(KBAGraphPartitionerTest, BricksHaveNoCyclicDependencies)
{
  // Non-uniform cuts divide a 7x5x4 grid of cells into 3x2x2 cellsets, which are partitioned one
  // per partition and then grouped onto fewer partitions. Each partition must be a brick of cells,
  // and no octant may have a cycle among partitions, as Koch, Baker and Alcouffe and Adams et al.
  // require.
  const auto grid = MakeCellGrid(7, 5, 4);
  // The analyzer reports an out-of-bounds read in boost::topological_sort's color map (a false
  // positive) and anchors it at this loop, the first user-code note on its path.
  // NOLINTNEXTLINE(clang-analyzer-security.ArrayBound)
  for (const int num_parts : {12, 6, 4, 3, 2, 1})
  {
    auto partitioner = MakeKBAWithCuts({1.0, 4.0}, {3.0}, {2.0});
    const auto pids = partitioner.Partition(grid.graph, grid.centroids, num_parts);
    EXPECT_EQ(std::set<int>(pids.begin(), pids.end()).size(), static_cast<size_t>(num_parts))
      << num_parts << " partitions";
    EXPECT_TRUE(PartitionsAreBricks(grid, pids, num_parts)) << num_parts << " partitions";
    EXPECT_TRUE(RankGraphAcyclicInEveryOctant(grid, pids, num_parts)) << num_parts << " partitions";
  }
}
