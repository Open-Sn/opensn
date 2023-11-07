#include "framework/physics/basic_options/basic_options.h"

const chi_physics::BasicOption&
chi_physics::BasicOptions::operator()(const std::string& option_name) const
{
  for (const auto& option : options_)
  {
    if (option.Name() == option_name) return option;
  }

  throw std::out_of_range("Basic option " + option_name + " does not appear to exist.");
}

const chi_physics::BasicOption&
chi_physics::BasicOptions::operator()(size_t index) const
{
  if (index < options_.size()) return options_[index];

  throw std::out_of_range("Basic option with index " + std::to_string(index) +
                          " does not appear to exist.");
}

chi_physics::BasicOption&
chi_physics::BasicOptions::operator[](const std::string& option_name)
{
  for (auto& option : options_)
  {
    if (option.Name() == option_name) return option;
  }

  throw std::out_of_range("Basic option \"" + option_name + "\" does not appear to exist.");
}

chi_physics::BasicOption&
chi_physics::BasicOptions::operator[](size_t index)
{
  if (index < options_.size()) return options_[index];

  throw std::out_of_range("Basic option with index " + std::to_string(index) +
                          " does not appear to exist.");
}

template <>
void
chi_physics::BasicOptions::AddOption<std::string>(const std::string& option_name,
                                                  const std::string& value)
{
  options_.emplace_back(option_name, value);
}

template <>
void
chi_physics::BasicOptions::AddOption<bool>(const std::string& option_name, const bool& value)
{
  options_.emplace_back(option_name, value);
}

template <>
void
chi_physics::BasicOptions::AddOption<int64_t>(const std::string& option_name, const int64_t& value)
{
  options_.emplace_back(option_name, value);
}

template <>
void
chi_physics::BasicOptions::AddOption<double>(const std::string& option_name, const double& value)
{
  options_.emplace_back(option_name, value);
}

size_t
chi_physics::BasicOptions::GetOptionIndexFromName(const std::string& option_name) const
{
  size_t index = 0;
  for (const auto& option : options_)
  {
    if (option.Name() == option_name) return index;
    ++index;
  }

  throw std::out_of_range("Basic option " + option_name + " does not appear to exist.");
}
