// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_base.h"
#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_base_mapping.h"

namespace opensn
{

/**
 * Generalization of the Galerkin Finite Element Method with piecewise linear basis functions for
 * use by a Continues Finite Element Method (CFEM).
 *
 * \ingroup doc_SpatialDiscretization
 */
class PieceWiseLinearContinuous : public PieceWiseLinearBase
{
public:
  /// Construct a shared object using the protected constructor.
  static std::shared_ptr<PieceWiseLinearContinuous>
  New(std::shared_ptr<MeshContinuum> grid, QuadratureOrder q_order = QuadratureOrder::SECOND);

  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            const UnknownManager& unknown_manager) const override;

  uint64_t MapDOF(const Cell& cell,
                  unsigned int node,
                  const UnknownManager& unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component) const override;

  uint64_t MapDOFLocal(const Cell& cell,
                       unsigned int node,
                       const UnknownManager& unknown_manager,
                       unsigned int unknown_id,
                       unsigned int component) const override;

  uint64_t MapDOF(const Cell& cell, unsigned int node) const override
  {
    return MapDOF(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }

  uint64_t MapDOFLocal(const Cell& cell, unsigned int node) const override
  {
    return MapDOFLocal(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }

  size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;

  std::vector<uint64_t> GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;

protected:
  /// Reorders the nodes for parallel computation in a Continuous Finite Element calculation.
  void OrderNodes();

  std::map<uint64_t, int64_t> node_mapping_;
  std::map<uint64_t, int64_t> ghost_node_mapping_;

private:
  explicit PieceWiseLinearContinuous(std::shared_ptr<MeshContinuum> grid, QuadratureOrder q_order);
};

} // namespace opensn
