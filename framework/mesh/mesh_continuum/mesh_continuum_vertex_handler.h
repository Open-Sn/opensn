// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <map>

namespace opensn
{

/// Manages a vertex map with custom calls.
class VertexHandler
{
  using GlobalIDMap = std::map<uint64_t, Vector3>;

public:
  // Iterators
  GlobalIDMap::iterator begin() { return m_global_id_vertex_map_.begin(); }
  GlobalIDMap::iterator end() { return m_global_id_vertex_map_.end(); }

  GlobalIDMap::const_iterator begin() const { return m_global_id_vertex_map_.begin(); }
  GlobalIDMap::const_iterator end() const { return m_global_id_vertex_map_.end(); }

  // Accessors
  Vector3& operator[](const uint64_t global_id) { return m_global_id_vertex_map_.at(global_id); }

  const Vector3& operator[](const uint64_t global_id) const
  {
    return m_global_id_vertex_map_.at(global_id);
  }

  // Utilities
  void Insert(const uint64_t global_id, const Vector3& vec)
  {
    m_global_id_vertex_map_.insert(std::make_pair(global_id, vec));
  }

  size_t GetNumLocallyStored() const { return m_global_id_vertex_map_.size(); }

  void Clear() { m_global_id_vertex_map_.clear(); }

private:
  std::map<uint64_t, Vector3> m_global_id_vertex_map_;
};

} // namespace opensn
