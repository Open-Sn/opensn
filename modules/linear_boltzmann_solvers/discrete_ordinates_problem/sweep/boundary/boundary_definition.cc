// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_definition.h"
#include <map>
#include <sstream>

namespace opensn
{

namespace
{

std::map<std::string, LBSBoundaryType> type_map = {{"vacuum", LBSBoundaryType::VACUUM},
                                                   {"isotropic", LBSBoundaryType::ISOTROPIC},
                                                   {"reflecting", LBSBoundaryType::REFLECTING},
                                                   {"arbitrary", LBSBoundaryType::ARBITRARY},
                                                   {"periodic", LBSBoundaryType::PERIODIC}};

std::string
UnsupportedParameterMessage(const std::string& boundary_name,
                            const std::string& boundary_type,
                            const std::string& param_name)
{
  std::ostringstream message;
  message << "Boundary '" << boundary_name << "' of type=\"" << boundary_type
          << "\" does not support \"" << param_name << "\".";
  return message.str();
}

std::string
MissingParameterMessage(const std::string& boundary_name,
                        const std::string& boundary_type,
                        const std::string& param_name)
{
  std::ostringstream message;
  message << "Boundary '" << boundary_name << "' of type=\"" << boundary_type << "\" requires \""
          << param_name << "\".";
  return message.str();
}

std::string
MissingEitherParameterMessage(const std::string& boundary_name,
                              const std::string& boundary_type,
                              const std::string& param_name1,
                              const std::string& param_name2)
{
  std::ostringstream message;
  message << "Boundary '" << boundary_name << "' of type=\"" << boundary_type
          << "\" requires either \"" << param_name1 << "\" or \"" << param_name2 << "\".";
  return message.str();
}

void
CheckForbiddenParams(const InputParameters& params,
                     const std::string& boundary_name,
                     const std::string& boundary_type,
                     const std::string& param_name)
{
  if (params.IsParameterValid(param_name))
    throw std::runtime_error(UnsupportedParameterMessage(boundary_name, boundary_type, param_name));
}

void
CheckRequiredParams(const InputParameters& params,
                    const std::string& boundary_name,
                    const std::string& boundary_type,
                    const std::string& param_name)
{
  if (not params.IsParameterValid(param_name))
    throw std::runtime_error(MissingParameterMessage(boundary_name, boundary_type, param_name));
}

void
CheckRequiredParams(const InputParameters& params,
                    const std::string& boundary_name,
                    const std::string& boundary_type,
                    const std::string& param_name1,
                    const std::string& param_name2)
{
  bool is_valid = (params.IsParameterValid(param_name1) xor params.IsParameterValid(param_name2));
  if (not is_valid)
    throw std::runtime_error(
      MissingEitherParameterMessage(boundary_name, boundary_type, param_name1, param_name2));
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
      CheckForbiddenParams(params, boundary_name, boundary_type, "time_function");
      CheckRequiredParams(params, boundary_name, boundary_type, "group_strength");
      params.RequireParameterBlockTypeIs("group_strength", ParameterBlockType::ARRAY);
      auto param_group_strength = params.GetParamVectorValue<double>("group_strength");
      if (param_group_strength.size() != num_groups)
      {
        std::ostringstream message;
        message << "Boundary '" << boundary_name
                << R"(' with type="isotropic" expected "group_strength" to contain )" << num_groups
                << " entries, one for each solver group, but received "
                << param_group_strength.size() << ".";
        throw std::runtime_error(message.str());
      }
      group_strength = std::move(param_group_strength);
      if (params.IsParameterValid("start_time"))
        start_time = params.GetParamValue<double>("start_time");
      if (params.IsParameterValid("end_time"))
        end_time = params.GetParamValue<double>("end_time");
      if (not(start_time <= end_time))
      {
        std::ostringstream message;
        message << "Boundary '" << boundary_name
                << "' with type=\"isotropic\" expected \"start_time\" to be less than or "
                   "equal to \"end_time\", but received start_time="
                << start_time << " and end_time=" << end_time << ".";
        throw std::runtime_error(message.str());
      }
      break;
    }
    case LBSBoundaryType::ARBITRARY:
    {
      CheckForbiddenParams(params, boundary_name, boundary_type, "group_strength");
      CheckForbiddenParams(params, boundary_name, boundary_type, "start_time");
      CheckForbiddenParams(params, boundary_name, boundary_type, "end_time");
      CheckRequiredParams(params, boundary_name, boundary_type, "function", "time_function");
      if (params.IsParameterValid("function"))
      {
        auto param_angular_flux_function =
          params.GetSharedPtrParam<AngularFluxFunction>("function", false);
        if (not param_angular_flux_function)
        {
          std::ostringstream message;
          message << "Boundary '" << boundary_name
                  << "' with type=\"arbitrary\" expected parameter \"function\" to hold a "
                     "non-null AngularFluxFunction.";
          throw std::runtime_error(message.str());
        }
        time_angular_flux_function = std::make_shared<AngularFluxTimeFunction>(
          [angular_flux_function =
             std::move(param_angular_flux_function)](int group, int direction, double time)
          { return (*angular_flux_function)(group, direction); });
      }
      else if (params.IsParameterValid("time_function"))
      {
        auto param_angular_flux_function =
          params.GetSharedPtrParam<AngularFluxTimeFunction>("time_function", false);
        if (not param_angular_flux_function)
        {
          std::ostringstream message;
          message << "Boundary '" << boundary_name
                  << "' with type=\"arbitrary\" expected parameter \"time_function\" to hold a "
                     "non-null AngularFluxFunction.";
          throw std::runtime_error(message.str());
        }
        time_angular_flux_function = std::move(param_angular_flux_function);
      }
      break;
    }
    default:
    {
      CheckForbiddenParams(params, boundary_name, boundary_type, "function");
      CheckForbiddenParams(params, boundary_name, boundary_type, "time_function");
      CheckForbiddenParams(params, boundary_name, boundary_type, "group_strength");
      CheckForbiddenParams(params, boundary_name, boundary_type, "start_time");
      CheckForbiddenParams(params, boundary_name, boundary_type, "end_time");
    }
  }
}

} // namespace opensn
