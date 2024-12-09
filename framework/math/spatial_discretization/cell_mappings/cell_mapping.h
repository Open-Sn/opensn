// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/vector.h"
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
  const Cell& ReferenceCell() const;

  /// Returns the grid on which the cell for this mapping lives.
  const MeshContinuum& ReferenceGrid() const;

  /// Returns the number of nodes on this element.
  size_t GetNumNodes() const;

  /// Returns the number of nodes on the given face.
  size_t GetNumFaceNodes(size_t face_index) const;

  const std::vector<std::vector<int>>& GetFaceNodeMappings() const;

  /// Returns the cell volume.
  double GetCellVolume() const;

  /// Returns the given face area.
  double GetFaceArea(size_t face_index) const;

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

  /// Returns the node locations associated with this element.
  const std::vector<Vector3>& GetNodeLocations() const;

  /// Makes the volumetric/internal finite element data for this element.
  virtual VolumetricFiniteElementData MakeVolumetricFiniteElementData() const = 0;

  /// Makes the surface finite element data for this element, at the specified face.
  virtual SurfaceFiniteElementData MakeSurfaceFiniteElementData(size_t face_index) const = 0;

  virtual ~CellMapping() = default;

protected:
  /**
   * This function gets called to compute the cell-volume andface-areas. If simple linear cells are
   * used then the default CellMapping::ComputeCellVolumeAndAreas can be used as this function.
   * Otherwise (i.e. for higher order elements, the child-class should bind a different function to
   * this.
   */
  using VolumeAndAreaFunction =
    std::function<void(const MeshContinuum&, const Cell&, double&, std::vector<double>&)>;

  CellMapping(const MeshContinuum& grid,
              const Cell& cell,
              size_t num_nodes,
              std::vector<Vector3> node_locations,
              std::vector<std::vector<int>> face_node_mappings,
              const VolumeAndAreaFunction& volume_area_function);

  /// Static method that all child elements can use as a default.
  static void ComputeCellVolumeAndAreas(const MeshContinuum& grid,
                                        const Cell& cell,
                                        double& volume,
                                        std::vector<double>& areas);

  const MeshContinuum& ref_grid_;
  const Cell& cell_;

  const size_t num_nodes_;
  const std::vector<Vector3> node_locations_;

  double volume_ = 0.0;
  std::vector<double> areas_;

  /**
   * For each cell face, map from the face node index to the corresponding cell node index. More
   * specifically, \p face_dof_mappings[f][fi], with \p fi the face node index of the face
   * identified by face index \p f, contains the corresponding cell node index.
   */
  const std::vector<std::vector<int>> face_node_mappings_;
};

} // namespace opensn
