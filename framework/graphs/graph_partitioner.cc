// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/graph_partitioner.h"

namespace opensn
{

InputParameters
GraphPartitioner::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  return params;
}

GraphPartitioner::GraphPartitioner(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
