// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

unsigned int
UnknownManager::AddUnknown(UnknownType unk_type, unsigned int dimension)
{
  unsigned int last_unknown_end = 0;
  if (not unknowns.empty())
    last_unknown_end = unknowns.back().GetMapEnd() + 1;

  if (unk_type == UnknownType::SCALAR)
  {
    unknowns.emplace_back(UnknownType::SCALAR, 1, last_unknown_end);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_2)
  {
    unknowns.emplace_back(UnknownType::VECTOR_2, 2, last_unknown_end);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_3)
  {
    unknowns.emplace_back(UnknownType::VECTOR_3, 3, last_unknown_end);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else if (unk_type == UnknownType::VECTOR_N)
  {
    if (dimension == 0)
    {
      std::ostringstream oss;
      oss << "UnknownManager: When adding unknown of type VECTOR_N, dimension cannot be 0";
      throw std::runtime_error(oss.str());
    }

    unknowns.emplace_back(UnknownType::VECTOR_N, dimension, last_unknown_end);
    unknowns.back().name = "Unknown_" + std::to_string(unknowns.size() - 1);
  }
  else
    throw std::runtime_error(
      "UnknownManager: Invalid call to AddUnknown. Unknown type not supported.");

  return unknowns.size();
}

unsigned int
UnknownManager::MapUnknown(int unknown_id, unsigned int component) const
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    std::ostringstream oss;
    oss << "UnknownManager: Call to MapUnknown failed (id = " << unknown_id << ")";
    throw std::runtime_error(oss.str());
  }
  return unknowns[unknown_id].GetMap(component);
}

unsigned int
UnknownManager::GetTotalUnknownStructureSize() const
{
  if (unknowns.empty())
    return 0;

  return unknowns.back().GetMapEnd() + 1;
}

void
UnknownManager::SetUnknownNumOffBlockConnections(int unknown_id, int num_conn)
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    std::ostringstream oss;
    oss << "UnknownManager: Call to SetUnknownNumOffBlockConnections failed (id = " << unknown_id
        << ")";
    throw std::runtime_error(oss.str());
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
    std::ostringstream oss;
    oss << "UnknownManager: Call to SetUnknownComponentNumOffBlockConnections failed "
        << "(id = " << unknown_id << ")";
    throw std::runtime_error(oss.str());
  }

  if (component < 0 or component >= unknowns[unknown_id].num_components)
  {
    std::ostringstream oss;
    oss << "UnknownManager: Call to SetUnknownComponentNumOffBlockConnections failed "
        << "(component id = " << component << ")";
    throw std::runtime_error(oss.str());
  }

  unknowns[unknown_id].num_off_block_connections[component] = num_conn;
}

void
UnknownManager::SetUnknownName(int unknown_id, const std::string& unk_name)
{
  if (unknown_id < 0 or unknown_id >= unknowns.size())
  {
    std::ostringstream oss;
    oss << "UnknownManager: Call to SetUnknownName failed (id = " << unknown_id << ")";
    throw std::runtime_error(oss.str());
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
    std::ostringstream oss;
    oss << "UnknownManager: Call to SetUnknownComponentName failed (id = " << unknown_id << ")";
    throw std::runtime_error(oss.str());
  }

  if (component < 0 or component >= unknowns[unknown_id].num_components)
  {
    std::ostringstream oss;
    oss << "UnknownManager: Call to SetUnknownComponentName failed "
        << "(component id = " << component << ")";
    throw std::runtime_error(oss.str());
  }

  unknowns[unknown_id].component_names[component] = component_name;
}

} // namespace opensn
