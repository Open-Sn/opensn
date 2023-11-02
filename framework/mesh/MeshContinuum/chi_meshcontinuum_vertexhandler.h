#ifndef CHI_MESHCONTINUUM_VERTEXHANDLER_H
#define CHI_MESHCONTINUUM_VERTEXHANDLER_H

#include "mesh/chi_meshvector.h"

#include <map>

namespace chi_mesh
{

/**Manages a vertex map with custom calls.*/
class VertexHandler
{
  typedef std::map<uint64_t, chi_mesh::Vector3> GlobalIDMap;

private:
  std::map<uint64_t, chi_mesh::Vector3> m_global_id_vertex_map;

public:
  // Iterators
  GlobalIDMap::iterator begin() { return m_global_id_vertex_map.begin(); }
  GlobalIDMap::iterator end() { return m_global_id_vertex_map.end(); }

  GlobalIDMap::const_iterator begin() const { return m_global_id_vertex_map.begin(); }
  GlobalIDMap::const_iterator end() const { return m_global_id_vertex_map.end(); }

  // Accessors
  chi_mesh::Vector3& operator[](const uint64_t global_id)
  {
    return m_global_id_vertex_map.at(global_id);
  }

  const chi_mesh::Vector3& operator[](const uint64_t global_id) const
  {
    return m_global_id_vertex_map.at(global_id);
  }

  // Utilities
  void Insert(const uint64_t global_id, const chi_mesh::Vector3& vec)
  {
    m_global_id_vertex_map.insert(std::make_pair(global_id, vec));
  }

  size_t NumLocallyStored() const { return m_global_id_vertex_map.size(); }

  void Clear() { m_global_id_vertex_map.clear(); }
};

} // namespace chi_mesh

#endif // CHI_MESHCONTINUUM_VERTEXHANDLER_H
