#pragma once

#ifdef OPENSN_WITH_LUA

#include <utility>

#include "framework/mesh/Cell/cell.h"

namespace lbs
{

struct ResponseFunctionDesignation
{
  const std::string name;
  const std::shared_ptr<chi_mesh::LogicalVolume> logical_volume;
  const std::string lua_functional;

  explicit ResponseFunctionDesignation(std::string in_name,
                                       std::shared_ptr<chi_mesh::LogicalVolume> in_logical_volume,
                                       std::string in_lua_function_name)
    : name(std::move(in_name)),
      logical_volume(std::move(in_logical_volume)),
      lua_functional(std::move(in_lua_function_name))
  {
  }

  std::vector<double> GetMGResponse(const chi_mesh::Cell& cell, size_t num_groups) const;
};

} // namespace lbs
#endif
