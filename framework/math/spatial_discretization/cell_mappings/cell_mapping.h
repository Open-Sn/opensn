// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector.h"
#include <memory>
#include <utility>
#include <vector>
#include <functional>

namespace opensn
{

class MeshContinuum;
struct Vector3;
class Cell;
class VolumetricFiniteElementData;
class SurfaceFiniteElementData;

/**
 * Base class for all cell mappings.
 *
 * \ingroup doc_CellMappings
 */
class CellMapping
{
public:
  /// Returns the cell this mapping is based on.
  const Cell& GetCell() const { return cell_; }

  /// Returns the grid on which the cell for this mapping lives.
  const std::shared_ptr<MeshContinuum> GetGrid() const { return grid_; }

  /// Returns the number of nodes on this element.
  size_t GetNumNodes() const { return num_nodes_; }

  /// Returns the number of nodes on the given face.
  size_t GetNumFaceNodes(size_t face_index) const;

  const std::vector<std::vector<int>>& GetFaceNodeMappings() const { return face_node_mappings_; }

  /// Returns the node locations associated with this element.
  const std::vector<Vector3>& GetNodeLocations() const { return node_locations_; }

  /**
   * Given the face index and the face node index, returns the index of the cell node the face node
   * corresponds to.
   */
  int MapFaceNode(size_t face_index, size_t face_node_index) const;

  /// Returns the value of the required shape function at the world xyz point.
  virtual double ShapeValue(int i, const Vector3& xyz) const = 0;

  /**
   * Populates all the shape function values at the given world xyz point. This method is optimized
   * to minimize reallocation of shape_values.
   */
  virtual void ShapeValues(const Vector3& xyz, Vector<double>& shape_values) const = 0;

  /// Returns the value of the required shape function gradient at the world xyz point.
  virtual Vector3 GradShapeValue(int i, const Vector3& xyz) const = 0;

  /**
   * Populates all the shape function gradient values at the given world xyz point. This method is
   * optimized to minimize reallocation of gradshape_values.
   */
  virtual void GradShapeValues(const Vector3& xyz,
                               std::vector<Vector3>& gradshape_values) const = 0;

  /// Makes the volumetric/internal finite element data for this element.
  virtual VolumetricFiniteElementData MakeVolumetricFiniteElementData() const = 0;

  /// Makes the surface finite element data for this element, at the specified face.
  virtual SurfaceFiniteElementData MakeSurfaceFiniteElementData(size_t face_index) const = 0;

  virtual ~CellMapping() = default;

protected:
  CellMapping(std::shared_ptr<MeshContinuum> grid,
              const Cell& cell,
              size_t num_nodes,
              std::vector<Vector3> node_locations,
              std::vector<std::vector<int>> face_node_mappings);

  const std::shared_ptr<MeshContinuum> grid_;
  const Cell& cell_;

  const size_t num_nodes_;
  const std::vector<Vector3> node_locations_;

  /**
   * For each cell face, map from the face node index to the corresponding cell node index. More
   * specifically, \p face_dof_mappings[f][fi], with \p fi the face node index of the face
   * identified by face index \p f, contains the corresponding cell node index.
   */
  const std::vector<std::vector<int>> face_node_mappings_;
};

} // namespace opensn
