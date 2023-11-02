#ifndef CHITECH_LAGRANGEDISCONTINUOUS_H
#define CHITECH_LAGRANGEDISCONTINUOUS_H

#include "LagrangeBase.h"

namespace chi_math::spatial_discretization
{

/**Generalization of the Galerkin Finite Element Method
 * with Lagrange basis functions
 * for use by a Discontinuous Finite Element Method (DFEM).
 * \ingroup doc_SpatialDiscretization*/
class LagrangeDiscontinuous : public LagrangeBase
{
public:
  // prevent anything else other than a shared pointer
  static std::shared_ptr<LagrangeDiscontinuous>
  New(const chi_mesh::MeshContinuum& grid,
      QuadratureOrder q_order = QuadratureOrder::SECOND,
      CoordinateSystemType cs_type = CoordinateSystemType::CARTESIAN);

  // 03
  void
  BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                       std::vector<int64_t>& nodal_nnz_off_diag,
                       const UnknownManager& unknown_manager) const override;

  // 04 Mappings
  int64_t MapDOF(const chi_mesh::Cell& cell,
                 unsigned int node,
                 const UnknownManager& unknown_manager,
                 unsigned int unknown_id,
                 unsigned int component) const override;

  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                      unsigned int node,
                      const UnknownManager& unknown_manager,
                      unsigned int unknown_id,
                      unsigned int component) const override;

  int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const override
  {
    return MapDOF(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }
  int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                      unsigned int node) const override
  {
    return MapDOFLocal(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }

  size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;

  std::vector<int64_t>
  GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;

protected:
  // 02
  void OrderNodes();

  std::vector<int64_t> cell_local_block_address_;
  std::vector<std::pair<uint64_t, int64_t>> neighbor_cell_block_address_;

private:
  // 00
  explicit LagrangeDiscontinuous(const chi_mesh::MeshContinuum& grid,
                                 QuadratureOrder q_order,
                                 CoordinateSystemType cs_type);
};

} // namespace chi_math::spatial_discretization

#endif // CHITECH_LAGRANGEDISCONTINUOUS_H
