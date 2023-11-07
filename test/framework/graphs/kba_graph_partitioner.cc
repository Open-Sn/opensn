#include "framework/console/console.h"

#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/object_factory.h"

#include "framework/mesh/mesh.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace chi_unit_tests
{

chi::ParameterBlock TestKBAGraphPartitioner00(const chi::InputParameters&);

RegisterWrapperFunction(chi_unit_tests,
                        TestKBAGraphPartitioner00,
                        nullptr,
                        TestKBAGraphPartitioner00);

chi::ParameterBlock
TestKBAGraphPartitioner00(const chi::InputParameters&)
{
  Chi::log.Log() << "GOLD_BEGIN";

  chi::ParameterBlock input_parameters;

  input_parameters.AddParameter("nx", 2);
  input_parameters.AddParameter("ny", 2);
  input_parameters.AddParameter("nz", 2);

  input_parameters.AddParameter("xcuts", std::vector<double>{0.0});
  input_parameters.AddParameter("ycuts", std::vector<double>{0.0});
  input_parameters.AddParameter("zcuts", std::vector<double>{0.0});

  chi::InputParameters valid_parameters = chi::KBAGraphPartitioner::GetInputParameters();

  valid_parameters.AssignParameters(input_parameters);

  chi::KBAGraphPartitioner partitioner(valid_parameters);

  std::vector<std::vector<uint64_t>> dummy_graph(8);
  std::vector<chi_mesh::Vector3> centroids = {{-1.0, -1.0, -1.0},
                                              {1.0, -1.0, -1.0},
                                              {-1.0, 1.0, -1.0},
                                              {1.0, 1.0, -1.0},
                                              {-1.0, -1.0, 1.0},
                                              {1.0, -1.0, 1.0},
                                              {-1.0, 1.0, 1.0},
                                              {1.0, 1.0, 1.0}};

  auto cell_pids = partitioner.Partition(dummy_graph, centroids, 2 * 2 * 2);

  for (const int64_t pid : cell_pids)
    Chi::log.Log() << pid;

  Chi::log.Log() << "GOLD_END";

  return chi::ParameterBlock();
}

} // namespace chi_unit_tests
