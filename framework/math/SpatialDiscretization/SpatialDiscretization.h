#pragma once

#include "opensn/framework/mesh/chi_mesh.h"
#include "opensn/framework/math/Quadratures/quadrature.h"
#include "opensn/framework/math/chi_math.h"
#include "opensn/framework/math/UnknownManager/unknown_manager.h"
#include "opensn/framework/mesh/Cell/cell.h"
#include "opensn/framework/math/SpatialDiscretization/CellMappings/CellMapping.h"

#include <petscksp.h>

#include <vector>
#include <map>

namespace chi_math
{
class SpatialDiscretization
{
public:
  const UnknownManager UNITARY_UNKNOWN_MANAGER;

  /**Utility method for getting node indices seperately for domain internal
   * local nodes, and boundary nodes.*/
  std::pair<std::set<uint32_t>, std::set<uint32_t>>
  MakeCellInternalAndBndryNodeIDs(const chi_mesh::Cell& cell) const;

  // 01 AddViewOfContinuum
  const CellMapping& GetCellMapping(const chi_mesh::Cell& cell) const;
  SpatialDiscretizationType Type() const;
  /**Returns the reference grid on which this discretization is based.*/
  const chi_mesh::MeshContinuum& Grid() const;
  CoordinateSystemType GetCoordinateSystemType() const;

  // 02 OrderNodes

  // 03
  /**Builds the sparsity pattern for a local block matrix compatible with
   * the given unknown manager. The modified vectors are: `nodal_nnz_in_diag`
   * which specifies for each row the number of non-zeros in the local diagonal
   * block, `nodal_nnz_off_diag` which specifies for each row the number of
   * non-zeros in the off-diagonal block.*/
  virtual void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                                    std::vector<int64_t>& nodal_nnz_off_diag,
                                    const UnknownManager& unknown_manager) const = 0;

  // 04 Mappings
  /**Maps the global address of a degree of freedom.*/
  virtual int64_t MapDOF(const chi_mesh::Cell& cell,
                         unsigned int node,
                         const UnknownManager& unknown_manager,
                         unsigned int unknown_id,
                         unsigned int component) const = 0;

  /**Maps the local address of a degree of freedom. This can include
   * ghost entries if the specific discretization has any.*/
  virtual int64_t MapDOFLocal(const chi_mesh::Cell& cell,
                              unsigned int node,
                              const UnknownManager& unknown_manager,
                              unsigned int unknown_id,
                              unsigned int component) const = 0;

  /**Maps the local address of a degree of freedom. This can include
   * ghost entries if the specific discretization has any. Default structure
   * here is a single scalar unknown.*/
  virtual int64_t MapDOF(const chi_mesh::Cell& cell, unsigned int node) const = 0;
  /**Maps the local address of a degree of freedom. This can include
   * ghost entries if the specific discretization has any. Default structure
   * here is a single scalar unknown.*/
  virtual int64_t MapDOFLocal(const chi_mesh::Cell& cell, unsigned int node) const = 0;

  // 05 Utils
  /**For the unknown structure in the unknown manager, returns the
   * number of local degrees-of-freedom.*/
  size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) const;
  /**For the unknown structure in the unknown manager, returns the
   * number of global degrees-of-freedom.*/
  size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) const;

  /**For the unknown structure in the unknown manager, returns the
   * number of ghost degrees-of-freedom.*/
  virtual size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const = 0;
  /**For the unknown structure in the unknown manager, returns the
   * global IDs of all the ghost degrees-of-freedom.*/
  virtual std::vector<int64_t> GetGhostDOFIndices(const UnknownManager& unknown_manager) const = 0;

  /**For the unknown structure in the unknown manager, returns the
   * number of local- and ghost degrees-of-freedom.*/
  size_t GetNumLocalAndGhostDOFs(const UnknownManager& unknown_manager) const;

  /**For the given cell, returns the number of relevant nodes. The same can
   * be achieved by retrieving the cell-to-element mapping first.*/
  size_t GetCellNumNodes(const chi_mesh::Cell& cell) const;

  /**For the given cell, returns a reference to the relevant node locations.
   * The same can be achieved by retrieving the cell-to-element mapping first.*/
  const std::vector<chi_mesh::Vector3>& GetCellNodeLocations(const chi_mesh::Cell& cell) const;

  /**\brief For each local cell, for each  face, for each face-node, provides a
  mapping to the adjacent cell's nodes.*/
  std::vector<std::vector<std::vector<int>>>
  MakeInternalFaceNodeMappings(double tolerance = 1.0e-12) const;

  /**Copies DOFs from the from_vector to the to_vector, however, both vectors
   * have the unknown-structure specified, and the given from_vec_uk_id is
   * copied to the to_vec_uk_id.*/
  void CopyVectorWithUnknownScope(const std::vector<double>& from_vector,
                                  std::vector<double>& to_vector,
                                  const UnknownManager& from_vec_uk_structure,
                                  unsigned int from_vec_uk_id,
                                  const UnknownManager& to_vec_uk_structure,
                                  unsigned int to_vec_uk_id) const;

  /**Develops a localized view of a petsc vector.
   * Each spatial discretization can have a specialization of this
   * method.*/
  virtual void LocalizePETScVector(Vec petsc_vector,
                                   std::vector<double>& local_vector,
                                   const UnknownManager& unknown_manager) const;
  /**Develops a localized view of a petsc vector.
   * Each spatial discretization can have a specialization of this
   * method.*/
  virtual void LocalizePETScVectorWithGhosts(Vec petsc_vector,
                                             std::vector<double>& local_vector,
                                             const UnknownManager& unknown_manager) const;

  /**Cartesian coordinate system spatial weighting function.*/
  static double CartesianSpatialWeightFunction(const chi_mesh::Vector3& point);
  /**Cylindrical coordinate system (RZ) spatial weighting function.*/
  static double CylindricalRZSpatialWeightFunction(const chi_mesh::Vector3& point);
  /**Spherical coordinate system (1D Spherical) spatial weighting function.*/
  static double Spherical1DSpatialWeightFunction(const chi_mesh::Vector3& point);

  typedef std::function<double(const chi_mesh::Vector3&)> SpatialWeightFunction;

  /**Returns the spatial weighting function appropriate to the discretization's
   * coordinate system.*/
  SpatialWeightFunction GetSpatialWeightingFunction() const;

  virtual ~SpatialDiscretization() = default;

protected:
  typedef SpatialDiscretizationType SDMType;
  // 00
  SpatialDiscretization(const chi_mesh::MeshContinuum& grid,
                        CoordinateSystemType cs_type,
                        SDMType sdm_type);

  const chi_mesh::MeshContinuum& ref_grid_;
  std::vector<std::unique_ptr<CellMapping>> cell_mappings_;
  std::map<uint64_t, std::shared_ptr<CellMapping>> nb_cell_mappings_;

  uint64_t local_block_address_ = 0;
  std::vector<uint64_t> locJ_block_address_;
  std::vector<uint64_t> locJ_block_size_;

  uint64_t local_base_block_size_ = 0;
  uint64_t globl_base_block_size_ = 0;

  const CoordinateSystemType coord_sys_type_;

private:
  const SpatialDiscretizationType type_;
};
} // namespace chi_math
