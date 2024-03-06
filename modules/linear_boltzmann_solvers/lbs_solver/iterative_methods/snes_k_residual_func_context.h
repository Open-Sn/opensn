#pragma once

#include <string>

namespace opensn
{
namespace lbs
{
struct KResidualFunctionContext
{
  std::string solver_name;
  double k_eff = 1.0;
};
} // namespace lbs
} // namespace opensn
