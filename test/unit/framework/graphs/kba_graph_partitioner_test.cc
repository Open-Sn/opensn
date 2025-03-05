#include "test/unit/opensn_unit_test.h"
#include "framework/graphs/kba_graph_partitioner.h"
#include <gmock/gmock.h>

using namespace opensn;

class KBAGraphPartitionerTest : public OpenSnUnitTest
{
};

TEST_F(KBAGraphPartitionerTest, Partition)
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
