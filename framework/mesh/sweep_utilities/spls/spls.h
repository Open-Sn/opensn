#pragma once

#include "framework/mesh/sweep_utilities/fluds/fluds.h"

namespace opensn
{

/**Contains the sweep plane data.*/
struct SPLS
{
  std::vector<int> item_id;
};

/**Stage-wise Task Dependency Graph.
 * Contains the global sweep plane data.*/
struct STDG
{
  std::vector<int> item_id;
};

} // namespace opensn
