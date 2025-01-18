#include "test/unit/opensn_unit_test.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"

void
OpenSnUnitTest::SetUp()
{
  opensn::Initialize();
}

void
OpenSnUnitTest::TearDown()
{
  opensn::Finalize();
}

std::shared_ptr<MeshContinuum>
OpenSnUnitTest::BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets) const
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
