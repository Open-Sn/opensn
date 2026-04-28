// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_definition.h"
#include <format>
#include <map>

namespace opensn
{

namespace
{

std::map<std::string, LBSBoundaryType> type_map = {{"vacuum", LBSBoundaryType::VACUUM},
                                                   {"isotropic", LBSBoundaryType::ISOTROPIC},
                                                   {"reflecting", LBSBoundaryType::REFLECTING},
                                                   {"arbitrary", LBSBoundaryType::ARBITRARY}};

void
CheckForbiddenParams(const InputParameters& params,
                     const std::string& boundary_name,
                     const std::string& boundary_type,
                     const std::string& param_name)
{
  if (params.IsParameterValid(param_name))
    throw std::runtime_error(std::format(R"(Boundary '{}' of type="{}" does not support "{}".)",
                                         boundary_name,
                                         boundary_type,
                                         param_name));
}

void
CheckRequiredParams(const InputParameters& params,
                    const std::string& boundary_name,
                    const std::string& boundary_type,
                    const std::string& param_name)
{
  if (not params.Has(param_name))
    throw std::runtime_error(std::format(
      R"(Boundary '{}' of type="{}" requires "{}".)", boundary_name, boundary_type, param_name));
}

} // namespace

BoundaryDefinition::BoundaryDefinition(const InputParameters& params, unsigned int num_groups)
{
  const auto boundary_name = params.GetParamValue<std::string>("name");
  const auto boundary_type = params.GetParamValue<std::string>("type");
  auto it = type_map.find(boundary_type);
  if (it == type_map.end())
  {
    throw std::runtime_error("Boundary \"" + boundary_name + "\" has unknown type \"" +
                             boundary_type + "\".");
  }
  type = it->second;

  switch (type)
  {
    case LBSBoundaryType::ISOTROPIC:
    {
      CheckForbiddenParams(params, boundary_name, boundary_type, "function");
      CheckRequiredParams(params, boundary_name, boundary_type, "group_strength");
      params.RequireParameterBlockTypeIs("group_strength", ParameterBlockType::ARRAY);
      auto param_group_strength = params.GetParamVectorValue<double>("group_strength");
      if (param_group_strength.size() != num_groups)
        throw std::runtime_error(std::format(
          "Boundary '{}' with type=\"isotropic\" expected \"group_strength\" to contain "
          "{} entries, one for each solver group, but received {}.",
          boundary_name,
          num_groups,
          param_group_strength.size()));
      group_strength = std::move(param_group_strength);
      break;
    }
    case LBSBoundaryType::ARBITRARY:
    {
      CheckForbiddenParams(params, boundary_name, boundary_type, "group_strength");
      CheckRequiredParams(params, boundary_name, boundary_type, "function");
      auto param_angular_flux_function =
        params.GetSharedPtrParam<AngularFluxFunction>("function", false);
      if (not param_angular_flux_function)
        throw std::runtime_error(std::format(
          "Boundary '{}' with type=\"arbitrary\" expected parameter \"function\" to hold a "
          "non-null AngularFluxFunction.",
          boundary_name));

      angular_flux_function = std::move(param_angular_flux_function);
      break;
    }
    default:
    {
      CheckForbiddenParams(params, boundary_name, boundary_type, "function");
      CheckForbiddenParams(params, boundary_name, boundary_type, "group_strength");
    }
  }
}

} // namespace opensn
