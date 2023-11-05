#include "opensn/framework/graphs/GraphPartitioner.h"

namespace chi
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

} // namespace chi
