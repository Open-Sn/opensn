#ifndef SPATIAL_DISCRETIZATION_FV_H
#define SPATIAL_DISCRETIZATION_FV_H

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/UnknownManager/unknown_manager.h"

#include <map>

//###################################################################
namespace chi_math::spatial_discretization
{
/**Spatial discretizations supporting Finite Volume representations.
\ingroup doc_SpatialDiscretization*/
class FiniteVolume : public SpatialDiscretization
{
private:
  std::map<uint64_t, uint64_t> neighbor_cell_local_ids_;

private:
  explicit FiniteVolume(const chi_mesh::MeshContinuum& grid, CoordinateSystemType cs_type);

public:
  virtual ~FiniteVolume() = default;
  // prevent anything else other than a shared pointer
  static std::shared_ptr<FiniteVolume>
  New(const chi_mesh::MeshContinuum& in_grid,
      CoordinateSystemType in_cs_type = CoordinateSystemType::CARTESIAN);

  // 01
  void CreateCellMappings();

  // 02 node ordering
protected:
  void OrderNodes();

  // 03 sparsity
public:
  void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                            std::vector<int64_t>& nodal_nnz_off_diag,
                            const UnknownManager& unknown_manager) const override;

  // 04a mappings
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
  int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const override
  {
    return MapDOFLocal(cell, node, UNITARY_UNKNOWN_MANAGER, 0, 0);
  }

  // 05 utils
  size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const override;
  std::vector<int64_t> GetGhostDOFIndices(const UnknownManager& unknown_manager) const override;
};

} // namespace chi_math::spatial_discretization

#endif