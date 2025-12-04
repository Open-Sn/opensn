// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <cstdint>

namespace opensn
{

class SPDS;

/**
 * \brief Get pointer to index data of a cell and the number of indexes for that cell.
 * \param device_node_indexes Pointer to the flatten node index structure on device.
 * \param cell_local_idx Cell local index.
 * \return Pointer to the indexes of the cell and the number of indexes (total number of face nodes)
 * for the cell.
 */
constexpr std::pair<const std::uint64_t*, std::uint64_t>
GetCellDataIndex(const std::uint64_t* device_node_indexes, const std::uint32_t& cell_local_idx)
{
  const std::uint64_t* cell_data =
    device_node_indexes + static_cast<std::uint64_t>(2 * cell_local_idx);
  return {device_node_indexes + cell_data[0], cell_data[1]};
}

struct FaceNodalMapping
{
  /// Face index on the neighbor cell.
  const int associated_face_;
  /// Face-node index on the neighbor face.
  const std::vector<short> face_node_mapping_;
  /// Cell-node index on the neighbor cell.
  const std::vector<short> cell_node_mapping_;

  FaceNodalMapping(int adj_face_idx,
                   const std::vector<short>& node_mapping,
                   const std::vector<short>& cell_node_mapping)
    : associated_face_(adj_face_idx),
      face_node_mapping_(node_mapping),
      cell_node_mapping_(cell_node_mapping)
  {
  }
};
using CellFaceNodalMapping = std::vector<FaceNodalMapping>;

class FLUDSCommonData
{
public:
  explicit FLUDSCommonData(const SPDS& spds,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  virtual ~FLUDSCommonData() = default;

  const SPDS& GetSPDS() const;
  const FaceNodalMapping& GetFaceNodalMapping(uint64_t cell_local_id, unsigned int face_id) const;

protected:
  const SPDS& spds_;
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings_;
};

} // namespace opensn
