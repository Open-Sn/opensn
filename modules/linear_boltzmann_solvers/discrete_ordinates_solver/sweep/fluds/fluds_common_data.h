// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <cstdint>

namespace opensn
{

class SPDS;

struct FaceNodalMapping
{
  /// Face index on the neighbor cell.
  const int associated_face_;
  /// Face-node index on the neighbor face.
  const std::vector<short> face_node_mapping_;
  /// Cell-node index on the neighbor cell.
  const std::vector<short> cell_node_mapping_;

  FaceNodalMapping(int ass_face,
                   const std::vector<short>& node_mapping,
                   const std::vector<short>& cell_node_mapping)
    : associated_face_(ass_face),
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
