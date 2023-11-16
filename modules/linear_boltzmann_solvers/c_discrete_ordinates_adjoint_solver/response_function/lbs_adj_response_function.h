#pragma once

#ifdef OPENSN_WITH_LUA

#include <utility>

#include "framework/mesh/cell/cell.h"

namespace opensn
{
namespace lbs
{

struct ResponseFunctionDesignation
{
  const std::string name;
  const std::shared_ptr<LogicalVolume> logical_volume;
  const std::string lua_functional;

  explicit ResponseFunctionDesignation(std::string in_name,
                                       std::shared_ptr<LogicalVolume> in_logical_volume,
                                       std::string in_lua_function_name)
    : name(std::move(in_name)),
      logical_volume(std::move(in_logical_volume)),
      lua_functional(std::move(in_lua_function_name))
  {
  }

  /** Calls the lua function associated with the response function and
   * returns a multigroup vector of the source values.*/
  std::vector<double> GetMGResponse(const Cell& cell, size_t num_groups) const;
};

} // namespace lbs
} // namespace opensn
#endif
