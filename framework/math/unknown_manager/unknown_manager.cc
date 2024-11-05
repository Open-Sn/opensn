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
  int last_unknown_end = -1;
  if (not unknowns.empty())
    last_unknown_end = unknowns.back().MapEnd();

  unsigned int new_unknown_index = unknowns.size();

  if (unk_type == UnknownType::SCALAR)
  {
    unknowns.emplace_back(UnknownType::SCALAR, 1, last_unknown_end + 1);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_2)
  {
    unknowns.emplace_back(UnknownType::VECTOR_2, 2, last_unknown_end + 1);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_3)
  {
    unknowns.emplace_back(UnknownType::VECTOR_3, 3, last_unknown_end + 1);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_N)
  {
    if (dimension == 0)
    {
      log.LogAllError() << "UnknownManager: When adding unknown of type VECTOR_N, "
                        << "the dimension must not be 0.";
      Exit(EXIT_FAILURE);
    }

    unknowns.emplace_back(UnknownType::VECTOR_N, dimension, last_unknown_end + 1);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
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
UnknownManager::MapUnknown(int unknown_id, unsigned int component) const
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.LogAllError() << "UnknownManager failed call to MapUnknown " << unknown_id;
    Exit(EXIT_FAILURE);
  }
  return unknowns[unknown_id].Map(component);
}

unsigned int
UnknownManager::TotalUnknownStructureSize() const
{
  if (unknowns.empty())
    return 0;

  return unknowns.back().MapEnd() + 1;
}

void
UnknownManager::SetUnknownNumOffBlockConnections(int unknown_id, int num_conn)
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownNumOffBlockConnections,"
                         " illegal index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  for (auto& val : unknowns[unknown_id].num_off_block_connections)
    val = num_conn;
}

void
UnknownManager::SetUnknownComponentNumOffBlockConnections(int unknown_id,
                                                          unsigned int component,
                                                          int num_conn)
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentNumOffBlockConnections,"
                         " illegal unknown index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  if (component < 0 or component >= unknowns[unknown_id].num_components)
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentNumOffBlockConnections,"
                         " illegal component index. "
                      << component;
    Exit(EXIT_FAILURE);
  }

  unknowns[unknown_id].num_off_block_connections[component] = num_conn;
}

void
UnknownManager::SetUnknownName(int unknown_id, const std::string& unk_name)
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownName,"
                         " illegal index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  unknowns[unknown_id].name = unk_name;
}

void
UnknownManager::SetUnknownComponentName(int unknown_id,
                                        unsigned int component,
                                        const std::string& component_name)
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentName,"
                         " illegal unknown index. "
                      << unknown_id;
    Exit(EXIT_FAILURE);
  }

  if (component < 0 or component >= unknowns[unknown_id].num_components)
  {
    log.LogAllError() << "UnknownManager failed call to SetUnknownComponentName,"
                         " illegal component index. "
                      << component;
    Exit(EXIT_FAILURE);
  }

  unknowns[unknown_id].component_names[component] = component_name;
}

} // namespace opensn
