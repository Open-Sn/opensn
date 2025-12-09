// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <string>
#include <stdexcept>

namespace opensn
{

/// Different types of variables.
enum class UnknownType
{
  SCALAR = 1,
  VECTOR_2 = 2,
  VECTOR_3 = 3,
  VECTOR_N = 4,
  TENSOR = 5
};

/// Nodal variable storage format.
enum class UnknownStorageType
{
  NODAL = 1,
  BLOCK = 2
};

/// Basic class for an variable.
class Unknown
{
public:
  explicit Unknown(UnknownType type, unsigned int num_components = 1, unsigned int map_begin = 0)
    : type(type),
      num_components((type == UnknownType::SCALAR)     ? 1
                     : (type == UnknownType::VECTOR_2) ? 2
                     : (type == UnknownType::VECTOR_3) ? 3
                                                       : num_components),
      map_begin(map_begin)
  {
    component_names.resize(num_components, std::string());
    for (unsigned int c = 0; c < num_components; ++c)
    {

      char buffer[100];
      snprintf(buffer, 100, " %03d", c);
      component_names[c] = buffer;
    }
    num_off_block_connections.resize(num_components, 0);
  }

  unsigned int GetMap(unsigned int component_number = 0) const
  {
    unsigned int map_value = 0;
    switch (type)
    {
      case UnknownType::SCALAR:
        if (component_number >= num_components)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=1"
                                  " for a SCALAR unknown.");
        map_value = 0;
        break;
      case UnknownType::VECTOR_2:
        if (component_number >= num_components)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=2"
                                  " for a VECTOR_2 unknown.");
        map_value = map_begin + component_number;
        break;
      case UnknownType::VECTOR_3:
        if (component_number >= num_components)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=3"
                                  " for a VECTOR_3 unknown.");
        map_value = map_begin + component_number;
        break;
      case UnknownType::VECTOR_N:
        if (component_number >= num_components)
          throw std::out_of_range(
            "Attempting to access component " + std::to_string(component_number) +
            ">=" + std::to_string(num_components) + " for a VECTOR_N unknown.");
        map_value = map_begin + component_number;
        break;
      case UnknownType::TENSOR:
        if (component_number >= num_components)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=" + std::to_string(num_components) + " for a TENSOR unknown.");
        map_value = map_begin + component_number;
        break;
      default:
        break;
    }

    return map_value;
  }
  unsigned int GetMapEnd() const { return map_begin + num_components - 1; }

  unsigned int GetNumComponents() const { return num_components; }

  const UnknownType type;
  const unsigned int num_components;
  const unsigned int map_begin;
  std::string name;
  std::vector<std::string> component_names;
  std::vector<int> num_off_block_connections;
};

/// General object for the management of unknowns in mesh-based mathematical model.
class UnknownManager
{
public:
  // Constructors
  explicit UnknownManager(UnknownStorageType storage_type = UnknownStorageType::NODAL) noexcept
    : dof_storage_type(storage_type)
  {
  }

  UnknownManager(std::initializer_list<std::pair<UnknownType, unsigned int>> unknown_info_list,
                 UnknownStorageType storage_type = UnknownStorageType::NODAL)
    : dof_storage_type(storage_type)
  {
    for (const auto& uk_info : unknown_info_list)
      AddUnknown(uk_info.first, uk_info.second);
  }

  explicit UnknownManager(const std::vector<Unknown>& unknown_info_list,
                          UnknownStorageType storage_type = UnknownStorageType::NODAL)
    : dof_storage_type(storage_type)
  {
    for (const auto& uk : unknown_info_list)
      AddUnknown(uk.type, uk.num_components);
  }

  UnknownManager(std::initializer_list<Unknown> unknowns,
                 UnknownStorageType storage_type = UnknownStorageType::NODAL)
    : dof_storage_type(storage_type)
  {
    std::size_t uk_id = 0;
    for (const auto& uk : unknowns)
    {
      AddUnknown(uk.type, uk.num_components);
      SetUnknownName(uk_id, uk.name);
      size_t comp_id = 0;
      for (const auto& comp_name : uk.component_names)
      {
        SetUnknownComponentName(uk_id, comp_id, comp_name);
        ++comp_id;
      }

      ++uk_id;
    }
  }

  UnknownManager(const UnknownManager& other) = default;
  UnknownManager& operator=(const UnknownManager& other) = default;

  std::size_t GetNumberOfUnknowns() const { return unknowns.size(); }
  const Unknown& GetUnknown(std::size_t id) const { return unknowns[id]; }

  void SetDOFStorageType(const UnknownStorageType storage_type) { dof_storage_type = storage_type; }

  UnknownStorageType GetDOFStorageType() const { return dof_storage_type; }

  void Clear() { unknowns.clear(); }

  /**
   * Adds an unknown to the manager. This method will figure out where the last unknown ends and
   * where to begin the next one.
   */
  std::size_t AddUnknown(UnknownType unk_type, unsigned int dimension = 0);

  /// Maps the unknown's component within the storage of a node.
  std::size_t MapUnknown(std::size_t unknown_id, unsigned int component = 0) const;

  /// Determines the total number of components over all unknowns.
  std::size_t GetTotalUnknownStructureSize() const;

  /**
   * Sets the number of off block connections for the given unknown. All the components will be set
   * to the same amount.
   */
  void SetUnknownNumOffBlockConnections(std::size_t unknown_id, int num_conn);

  /// Sets the number of off block connections for the given unknown-component pair.
  void SetUnknownComponentNumOffBlockConnections(std::size_t unknown_id,
                                                 unsigned int component,
                                                 int num_conn);

  /// Sets a unk_name for the indicated unknown
  void SetUnknownName(std::size_t unknown_id, const std::string& unk_name);

  /// Sets the component_name to be associated with each component of the unknown.
  void SetUnknownComponentName(std::size_t unknown_id,
                               unsigned int component,
                               const std::string& component_name);

  ~UnknownManager() = default;

  std::vector<Unknown> unknowns;
  UnknownStorageType dof_storage_type;

public:
  // Utilities
  static UnknownManager GetUnitaryUnknownManager()
  {
    return UnknownManager({std::make_pair(UnknownType::SCALAR, 0)});
  }
};

} // namespace opensn
