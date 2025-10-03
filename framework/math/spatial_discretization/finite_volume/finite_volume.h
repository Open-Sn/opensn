// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include <map>

namespace opensn
{

/**
 * Spatial discretizations supporting Finite Volume representations.
 *
 * \ingroup doc_SpatialDiscretization
 */
class FiniteVolume : public SpatialDiscretization
{
private:
  std::map<uint64_t, uint64_t> neighbor_cell_local_ids_;

private:
  explicit FiniteVolume(std::shared_ptr<MeshContinuum> grid);

public:
  ~FiniteVolume() override = default;

  /// Publicly accessible construction handler.
  static std::shared_ptr<FiniteVolume> New(std::shared_ptr<MeshContinuum> grid);

  void CreateCellMappings();

protected:
  /// Develops node ordering per location.
  void OrderNodes();

public:
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            const UnknownManager& unknown_manager) const override;

  uint64_t MapDOF(const Cell& cell,
                  unsigned int node,
                  const UnknownManager& unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component) const override;

  int64_t MapDOFLocal(const Cell& cell,
                      unsigned int node,
                      const UnknownManager& unknown_manager,
                      unsigned int unknown_id,
                      unsigned int component) const override;

  uint64_t MapDOF(const Cell& cell, unsigned int node) const override
  {
    return MapDOF(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }

  int64_t MapDOFLocal(const Cell& cell, unsigned int node) const override
  {
    return MapDOFLocal(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }

  size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;
  std::vector<uint64_t> GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;
};

} // namespace opensn
