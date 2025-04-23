// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <petscksp.h>
#include <vector>
#include <map>
#include <set>

namespace opensn
{
class SpatialDiscretization
{
public:
  const UnknownManager UNITARY_UNKNOWN_MANAGER;

  /**
   * Utility method for getting node indices seperately for domain internal local nodes, and
   * boundary nodes.
   */
  std::pair<std::set<uint32_t>, std::set<uint32_t>>
  MakeCellInternalAndBndryNodeIDs(const Cell& cell) const;

  const CellMapping& GetCellMapping(const Cell& cell) const;
  SpatialDiscretizationType GetType() const;

  /// Returns the reference grid on which this discretization is based.
  const std::shared_ptr<MeshContinuum> GetGrid() const;

  /**
   * Builds the sparsity pattern for a local block matrix compatible withthe given unknown manager.
   * The modified vectors are: `nodal_nnz_in_diag` which specifies for each row the number of
   * non-zeros in the local diagonal block, `nodal_nnz_off_diag` which specifies for each row the
   * number of non-zeros in the off-diagonal block.
   */
  virtual void BuildSparsityPattern(std::vector<int64_t>& nodal_nnz_in_diag,
                                    std::vector<int64_t>& nodal_nnz_off_diag,
                                    const UnknownManager& unknown_manager) const = 0;

  /// Maps the global address of a degree of freedom.
  virtual int64_t MapDOF(const Cell& cell,
                         unsigned int node,
                         const UnknownManager& unknown_manager,
                         unsigned int unknown_id,
                         unsigned int component) const = 0;

  /// Maps the local address of a degree of freedom. This can include ghost entries if the specific
  /// discretization has any.
  virtual int64_t MapDOFLocal(const Cell& cell,
                              unsigned int node,
                              const UnknownManager& unknown_manager,
                              unsigned int unknown_id,
                              unsigned int component) const = 0;

  /**
   * Maps the local address of a degree of freedom. This can include ghost entries if the specific
   * discretization has any. Default structure here is a single scalar unknown.
   */
  virtual int64_t MapDOF(const Cell& cell, unsigned int node) const = 0;

  /**
   * Maps the local address of a degree of freedom. This can include ghost entries if the specific
   * discretization has any. Default structure here is a single scalar unknown.
   */
  virtual int64_t MapDOFLocal(const Cell& cell, unsigned int node) const = 0;

  /// Returns the number of local nodes used in this discretization.
  size_t GetNumLocalNodes() const;

  /// Returns the number of global nodes used in this discretization.
  size_t GetNumGlobalNodes() const;

  /**
   * For the unknown structure in the unknown manager, returns the number of local
   * degrees-of-freedom.
   */
  size_t GetNumLocalDOFs(const UnknownManager& unknown_manager) const;

  /**
   * For the unknown structure in the unknown manager, returns the number of global
   * degrees-of-freedom.
   */
  size_t GetNumGlobalDOFs(const UnknownManager& unknown_manager) const;

  /**
   * For the unknown structure in the unknown manager, returns the number of ghost
   * degrees-of-freedom.
   */
  virtual size_t GetNumGhostDOFs(const UnknownManager& unknown_manager) const = 0;

  /**
   *For the unknown structure in the unknown manager, returns the global IDs of all the ghost
   * degrees-of-freedom.
   */
  virtual std::vector<int64_t> GetGhostDOFIndices(const UnknownManager& unknown_manager) const = 0;

  /**
   * For the unknown structure in the unknown manager, returns the number of local- and ghost
   * degrees-of-freedom.
   */
  size_t GetNumLocalAndGhostDOFs(const UnknownManager& unknown_manager) const;

  /**
   * For the given cell, returns the number of relevant nodes. The same can be achieved by
   * retrieving the cell-to-element mapping first.
   */
  size_t GetCellNumNodes(const Cell& cell) const;

  /**
   * For the given cell, returns a reference to the relevant node locations. The same can be
   * achieved by retrieving the cell-to-element mapping first.
   */
  const std::vector<Vector3>& GetCellNodeLocations(const Cell& cell) const;

  /**
   * For each cell, for each face of that cell, for each node on that face, maps to which local
   * node on the adjacent cell that node position corresponds.
   *
   * \param tolerance double. Tolerance to use to determine if two node locations are
   *                            equal. [Default: 1.0e-12]
   *
   *  For example consider two adjacent quadrilaterals.
   *
   * \verbatim
   * o--------o--------o                       o--------o
   * |3      2|3      2|                       |    2   |
   * |  101   |   102  | , face ids for both:  |3      1|
   * |0      1|0      1|                       |    0   |
   * o--------o--------o                       o--------o
   * internal face for cell 101 is face-1, ccw orientated
   * --o
   *  1|
   *   |
   *  0|
   * --o
   * internal face for cell 102 is face-3, ccw orientated
   * o-
   * |0
   * |
   * |1
   * o-
   *
   * mapping[101][1][0] = 0
   * mapping[101][1][1] = 3
   *
   * mapping[102][3][0] = 2
   * mapping[102][3][1] = 1
   * \endverbatim
   */
  std::vector<std::vector<std::vector<int>>>
  MakeInternalFaceNodeMappings(double tolerance = 1.0e-12) const;

  /**
   * Copy part of vector A to vector B. Suppose vector A's entries are managed `UnknownManager` A
   * (`uk_manA`) and that the entries of the vector B are managed by `UnknownManager` B (`uk_manB`).
   * This function copies the entries associated with an unknown with id `uk_id_A` in `uk_manA` from
   * vector A to vector B such that the entries in vector B are aligned with the entries of an
   * unknown with id `uk_id_B` in `uk_manB`. All the components are copied.
   *
   * \param from_vector Vector to copy from.
   * \param to_vector Vector to copy to.
   * \param from_vec_uk_structure Unknown manager for vector A.
   * \param from_vec_uk_id Unknown-id in unknown manager A.
   * \param to_vec_uk_structure Unknown manager for vector B.
   * \param to_vec_uk_id Unknown-id in unknown manager B.
   */
  void CopyVectorWithUnknownScope(const std::vector<double>& from_vector,
                                  std::vector<double>& to_vector,
                                  const UnknownManager& from_vec_uk_structure,
                                  unsigned int from_vec_uk_id,
                                  const UnknownManager& to_vec_uk_structure,
                                  unsigned int to_vec_uk_id) const;

  /**
   * Develops a localized view of a petsc vector. Each spatial discretization can have a
   * specialization of this method.
   */
  virtual void LocalizePETScVector(Vec petsc_vector,
                                   std::vector<double>& local_vector,
                                   const UnknownManager& unknown_manager) const;
  /**
   * Develops a localized view of a petsc vector. Each spatial discretization can have a
   * specialization of this method.
   */
  virtual void LocalizePETScVectorWithGhosts(Vec petsc_vector,
                                             std::vector<double>& local_vector,
                                             const UnknownManager& unknown_manager) const;

  /// Cartesian coordinate system spatial weighting function.
  static double CartesianSpatialWeightFunction(const Vector3& point);

  /// Cylindrical coordinate system (RZ) spatial weighting function.
  static double CylindricalRZSpatialWeightFunction(const Vector3& point);

  /// Spherical coordinate system (1D Spherical) spatial weighting function.
  static double Spherical1DSpatialWeightFunction(const Vector3& point);

  using SpatialWeightFunction = std::function<double(const Vector3&)>;

  /// Returns the spatial weighting function appropriate to the discretization's coordinate system.
  SpatialWeightFunction GetSpatialWeightingFunction() const;

  virtual ~SpatialDiscretization() = default;

protected:
  SpatialDiscretization(const std::shared_ptr<MeshContinuum> grid,
                        SpatialDiscretizationType sdm_type);

  const std::shared_ptr<MeshContinuum> grid_;
  std::vector<std::unique_ptr<CellMapping>> cell_mappings_;
  std::map<uint64_t, std::shared_ptr<CellMapping>> nb_cell_mappings_;

  uint64_t local_block_address_ = 0;
  std::vector<uint64_t> locJ_block_address_;
  std::vector<uint64_t> locJ_block_size_;

  uint64_t local_base_block_size_ = 0;
  uint64_t global_base_block_size_ = 0;

private:
  const SpatialDiscretizationType type_;
};

} // namespace opensn
