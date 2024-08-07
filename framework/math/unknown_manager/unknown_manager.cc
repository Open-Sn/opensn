// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

unsigned int
UnknownManager::AddUnknown(UnknownType unk_type, unsigned int dimension)
{
  auto& log = Logger::GetInstance();

  unsigned int last_unknown_end = -1;
  if (not unknowns_.empty())
    last_unknown_end = unknowns_.back().GetMapEnd();

  unsigned int new_unknown_index = unknowns_.size();

  if (unk_type == UnknownType::SCALAR)
  {
    unknowns_.emplace_back(UnknownType::SCALAR, 1, last_unknown_end + 1);
    unknowns_.back().name_ = "Unknown_" + std::to_string(unknowns_.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_2)
  {
    unknowns_.emplace_back(UnknownType::VECTOR_2, 2, last_unknown_end + 1);
    unknowns_.back().name_ = "Unknown_" + std::to_string(unknowns_.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_3)
  {
    unknowns_.emplace_back(UnknownType::VECTOR_3, 3, last_unknown_end + 1);
    unknowns_.back().name_ = "Unknown_" + std::to_string(unknowns_.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_N)
  {
    if (dimension == 0)
    {
      log.LogAllError() << "UnknownManager: When adding unknown of type VECTOR_N, "
                        << "the dimension must not be 0.";
      Exit(EXIT_FAILURE);
    }

    unknowns_.emplace_back(UnknownType::VECTOR_N, dimension, last_unknown_end + 1);
    unknowns_.back().name_ = "Unknown_" + std::to_string(unknowns_.size() - 1);
  }
  else if (unk_type == UnknownType::TENSOR)
  {
    if (dimension == 0 or dimension == 1)
    {
      log.LogAllError() << "UnknownManager: When adding unknown of type TENSOR, "
                        << "the dimension must not be 0 or 1.";
      Exit(EXIT_FAILURE);
    }

    throw std::invalid_argument("UnknownManager: TENSOR unknowns are not "
                                "supported yet.");
  }
  else
  {
    throw std::logic_error("UnknownManager: Invalid call to AddUnknown(). "
                           "Unknown type is probably not supported yet.");
  }

  return new_unknown_index;
}

unsigned int
UnknownManager::MapUnknown(unsigned int unknown_id, unsigned int component) const
{
  auto& log = Logger::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns_.size())
  {
    log.LogAllError() << "UnknownManager failed call to MapUnknown " << unknown_id;
    Exit(EXIT_FAILURE);
  }
  return unknowns_[unknown_id].GetMap(component);
}

unsigned int
UnknownManager::GetTotalUnknownStructureSize() const
{
  if (unknowns_.empty())
    return 0;

  return unknowns_.back().GetMapEnd() + 1;
}

void
UnknownManager::SetUnknownNumOffBlockConnections(unsigned int unknown_id, int num_conn)
{
  auto& log = Logger::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns_.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownNumOffBlockConnections,"
                         " illegal index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  for (auto& val : unknowns_[unknown_id].num_off_block_connections_)
    val = num_conn;
}

void
UnknownManager::SetUnknownComponentNumOffBlockConnections(unsigned int unknown_id,
                                                          unsigned int component,
                                                          int num_conn)
{
  auto& log = Logger::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns_.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentNumOffBlockConnections,"
                         " illegal unknown index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  if (component < 0 or component >= unknowns_[unknown_id].num_components_)
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentNumOffBlockConnections,"
                         " illegal component index. "
                      << component;
    Exit(EXIT_FAILURE);
  }

  unknowns_[unknown_id].num_off_block_connections_[component] = num_conn;
}

void
UnknownManager::SetUnknownTextName(unsigned int unknown_id, const std::string& name)
{
  auto& log = Logger::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns_.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownTextName,"
                         " illegal index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  unknowns_[unknown_id].name_ = name;
}

void
UnknownManager::SetUnknownComponentName(unsigned int unknown_id,
                                        unsigned int component,
                                        const std::string& name)
{
  auto& log = Logger::GetInstance();

  if (unknown_id < 0 or unknown_id >= unknowns_.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentName,"
                         " illegal unknown index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  if (component < 0 or component >= unknowns_[unknown_id].num_components_)
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentName,"
                         " illegal component index. "
                      << component;
    Exit(EXIT_FAILURE);
  }

  unknowns_[unknown_id].component_names_[component] = name;
}

} // namespace opensn
