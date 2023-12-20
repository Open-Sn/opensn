#include "framework/graphs/graph_partitioner.h"

namespace opensn
{

InputParameters
GraphPartitioner::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  return params;
}

GraphPartitioner::GraphPartitioner(const InputParameters& params) : ChiObject(params)
{
}

} // namespace opensn
