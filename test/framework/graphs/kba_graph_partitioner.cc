#include "lua/framework/console/console.h"
#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/object_factory.h"

#include "framework/mesh/mesh.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock TestKBAGraphPartitioner00(const InputParameters&);

RegisterWrapperFunctionNamespace(unit_tests,
                                 TestKBAGraphPartitioner00,
                                 nullptr,
                                 TestKBAGraphPartitioner00);

ParameterBlock
TestKBAGraphPartitioner00(const InputParameters&)
{
  opensn::log.Log() << "GOLD_BEGIN";

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

  for (const int64_t pid : cell_pids)
    opensn::log.Log() << pid;

  opensn::log.Log() << "GOLD_END";

  return ParameterBlock();
}

} //  namespace unit_tests
