#pragma once

#include <utility>

#include "framework/mesh/cell/cell.h"

namespace opensn
{
class ResponseFunction;

namespace lbs
{

struct ResponseFunctionDesignation
{
  const std::string name;
  const std::shared_ptr<LogicalVolume> logical_volume;
  std::shared_ptr<ResponseFunction> response_function;

  explicit ResponseFunctionDesignation(std::string in_name,
                                       std::shared_ptr<LogicalVolume> in_logical_volume,
                                       std::shared_ptr<ResponseFunction> in_response_function)
    : name(std::move(in_name)),
      logical_volume(std::move(in_logical_volume)),
      response_function(in_response_function)
  {
  }

  /** Calls the lua function associated with the response function and
   * returns a multigroup vector of the source values.*/
  std::vector<double> GetMGResponse(const Cell& cell, size_t num_groups) const;
};

} // namespace lbs
} // namespace opensn
