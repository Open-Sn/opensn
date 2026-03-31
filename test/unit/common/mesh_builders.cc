#include "test/unit/common/mesh_builders.h"
#include "framework/parameters/parameter_block.h"
#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"

using namespace opensn;

std::shared_ptr<MeshContinuum>
BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets)
{
  ParameterBlock array("node_sets");
  for (std::size_t i = 0; i < node_sets.size(); ++i)
    array.AddParameter(ParameterBlock(std::to_string(i + 1), node_sets[i]));
  array.ChangeToArray();

  ParameterBlock block("");
  block.AddParameter(array);

  auto params = OrthogonalMeshGenerator::GetInputParameters();
  params.AssignParameters(block);
  return OrthogonalMeshGenerator(params).Execute();
}
