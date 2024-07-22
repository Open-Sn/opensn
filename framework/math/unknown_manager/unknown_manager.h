// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <string>
#include <stdexcept>

namespace opensn
{

/**Different types of variables.*/
enum class UnknownType
{
  SCALAR = 1,
  VECTOR_2 = 2,
  VECTOR_3 = 3,
  VECTOR_N = 4,
  TENSOR = 5
};

/**Nodal variable storage format.*/
enum class UnknownStorageType
{
  NODAL = 1,
  BLOCK = 2
};

/**Basic class for an variable.*/
class Unknown
{
public:
  const UnknownType type_;
  const unsigned int num_components_;
  const unsigned int map_begin_;
  std::string text_name_;
  std::vector<std::string> component_text_names_;
  std::vector<int> num_off_block_connections_;

public:
  explicit Unknown(UnknownType type, unsigned int num_components = 1, unsigned int map_begin = 0)
    : type_(type),
      num_components_((type_ == UnknownType::SCALAR)     ? 1
                      : (type_ == UnknownType::VECTOR_2) ? 2
                      : (type_ == UnknownType::VECTOR_3) ? 3
                                                         : num_components),
      map_begin_(map_begin)
  {
    component_text_names_.resize(num_components_, std::string());
    for (unsigned int c = 0; c < num_components_; ++c)
    {

      char buffer[100];
      snprintf(buffer, 100, " %03d", c);
      component_text_names_[c] = buffer;
    }
    num_off_block_connections_.resize(num_components_, 0);
  }

  unsigned int GetMap(unsigned int component_number = 0) const
  {
    unsigned int map_value = 0;
    switch (type_)
    {
      case UnknownType::SCALAR:
        if (component_number >= num_components_)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=1"
                                  " for a SCALAR unknown.");
        map_value = 0;
        break;
      case UnknownType::VECTOR_2:
        if (component_number >= num_components_)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=2"
                                  " for a VECTOR_2 unknown.");
        map_value = map_begin_ + component_number;
        break;
      case UnknownType::VECTOR_3:
        if (component_number >= num_components_)
          throw std::out_of_range("Attempting to access component " +
                                  std::to_string(component_number) +
                                  ">=3"
                                  " for a VECTOR_3 unknown.");
        map_value = map_begin_ + component_number;
        break;
      case UnknownType::VECTOR_N:
        if (component_number >= num_components_)
          throw std::out_of_range(
            "Attempting to access component " + std::to_string(component_number) +
            ">=" + std::to_string(num_components_) + " for a VECTOR_N unknown.");
        map_value = map_begin_ + component_number;
        break;
      case UnknownType::TENSOR:
        if (component_number >= num_components_)
          throw std::out_of_range(
            "Attempting to access component " + std::to_string(component_number) +
            ">=" + std::to_string(num_components_) + " for a TENSOR unknown.");
        map_value = map_begin_ + component_number;
        break;
      default:
        break;
    }

    return map_value;
  }
  unsigned int GetMapEnd() const { return map_begin_ + num_components_ - 1; }

  unsigned int NumComponents() const { return num_components_; }
};

/**General object for the management of unknowns in mesh-based
 * mathematical model.*/
class UnknownManager
{
private:
public:
  std::vector<Unknown> unknowns_;
  UnknownStorageType dof_storage_type_;

  // Constructors
  explicit UnknownManager(UnknownStorageType storage_type = UnknownStorageType::NODAL) noexcept
    : dof_storage_type_(storage_type)
  {
  }

  UnknownManager(std::initializer_list<std::pair<UnknownType, unsigned int>> unknown_info_list,
                 UnknownStorageType storage_type = UnknownStorageType::NODAL) noexcept
    : dof_storage_type_(storage_type)
  {
    for (const auto& uk_info : unknown_info_list)
      AddUnknown(uk_info.first, uk_info.second);
  }

  explicit UnknownManager(const std::vector<Unknown>& unknown_info_list,
                          UnknownStorageType storage_type = UnknownStorageType::NODAL) noexcept
    : dof_storage_type_(storage_type)
  {
    for (const auto& uk : unknown_info_list)
      AddUnknown(uk.type_, uk.num_components_);
  }

  UnknownManager(std::initializer_list<Unknown> unknowns,
                 UnknownStorageType storage_type = UnknownStorageType::NODAL) noexcept
    : dof_storage_type_(storage_type)
  {
    size_t ukid = 0;
    for (const auto& uk : unknowns)
    {
      AddUnknown(uk.type_, uk.num_components_);
      SetUnknownTextName(ukid, uk.text_name_);
      size_t compid = 0;
      for (const auto& comp_text_name : uk.component_text_names_)
      {
        SetUnknownComponentTextName(ukid, compid, comp_text_name);
        ++compid;
      }

      ++ukid;
    }
  }

  UnknownManager(const UnknownManager& other) = default;
  UnknownManager& operator=(const UnknownManager& other) = default;

  // Utilities
  static UnknownManager GetUnitaryUnknownManager()
  {
    return UnknownManager({std::make_pair(UnknownType::SCALAR, 0)});
  }

  size_t NumberOfUnknowns() const { return unknowns_.size(); }
  const Unknown& GetUnknown(size_t id) const { return unknowns_[id]; }

  void SetDOFStorageType(const UnknownStorageType storage_type)
  {
    dof_storage_type_ = storage_type;
  }

  UnknownStorageType GetDOFStorageType() const { return dof_storage_type_; }

  void Clear() { unknowns_.clear(); }

  /**Adds an unknown to the manager. This method will figure out
   * where the last unknown ends and where to begin the next one.*/
  unsigned int AddUnknown(UnknownType unk_type, unsigned int dimension = 0);

  /**Maps the unknown's component within the storage of a node.*/
  unsigned int MapUnknown(unsigned int unknown_id, unsigned int component = 0) const;

  /**Determines the total number of components over all unknowns.*/
  unsigned int GetTotalUnknownStructureSize() const;

  /**Sets the number of off block connections for the given unknown.
   * All the components will be set to the same amount.*/
  void SetUnknownNumOffBlockConnections(unsigned int unknown_id, int num_conn);

  /**Sets the number of off block connections for the given unknown-
   * component pair.*/
  void SetUnknownComponentNumOffBlockConnections(unsigned int unknown_id,
                                                 unsigned int component,
                                                 int num_conn);

  /**Sets a text name for the indicated unknown.*/
  void SetUnknownTextName(unsigned int unknown_id, const std::string& text_name);

  /**Sets the text name to be associated with each component of the
   * unknown.*/
  void SetUnknownComponentTextName(unsigned int unknown_id,
                                   unsigned int component,
                                   const std::string& text_name);

  ~UnknownManager() = default;
};

} // namespace opensn
