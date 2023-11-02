#ifndef CHITECH_CELLMAPPING_H
#define CHITECH_CELLMAPPING_H

#include <memory>
#include <utility>
#include <vector>
#include <functional>

// ################################################################### Fwd
// Decls.
namespace chi_mesh
{
class MeshContinuum;
struct Vector3;
class Cell;
} // namespace chi_mesh

namespace chi_math::finite_element
{
class VolumetricQuadraturePointData;
class SurfaceQuadraturePointData;
} // namespace chi_math::finite_element

namespace chi_math
{
// ################################################################### Class def
/**Base class for all cell mappings.
* \ingroup doc_CellMappings*/
class CellMapping
{
public:
  // 00
  /**Returns the cell this mapping is based on.*/
  const chi_mesh::Cell& ReferenceCell() const;

  /**Returns the grid on which the cell for this mapping lives.*/
  const chi_mesh::MeshContinuum& ReferenceGrid() const;

  /**Returns the number of nodes on this element.*/
  size_t NumNodes() const;
  /**Returns the number of nodes on the given face.*/
  size_t NumFaceNodes(size_t face_index) const;

  const std::vector<std::vector<int>>& GetFaceNodeMappings() const;

  /**Returns the cell volume.*/
  double CellVolume() const;
  /**Returns the given face area.*/
  double FaceArea(size_t face_index) const;

  /**Given the face index and the face node index, returns the index
   * of the cell node the face node corresponds to.*/
  int MapFaceNode(size_t face_index, size_t face_node_index) const;

  // 02 ShapeFuncs
  /**Returns the value of the required shape function at the world xyz point.*/
  virtual double ShapeValue(int i, const chi_mesh::Vector3& xyz) const = 0;

  /**Populates all the shape function values at the given world xyz point. This
  * method is optimized to minimize reallocation of shape_values.*/
  virtual void ShapeValues(const chi_mesh::Vector3& xyz,
                           std::vector<double>& shape_values) const = 0;

  /**Returns the value of the required shape function gradient at the world xyz
   * point.*/
  virtual chi_mesh::Vector3
  GradShapeValue(int i, const chi_mesh::Vector3& xyz) const = 0;

  /**Populates all the shape function gradient values at the given world xyz
  * point. This method is optimized to minimize reallocation of
  * gradshape_values.*/
  virtual void
  GradShapeValues(const chi_mesh::Vector3& xyz,
                  std::vector<chi_mesh::Vector3>& gradshape_values) const = 0;

  /**Returns the node locations associated with this element.*/
  const std::vector<chi_mesh::Vector3>& GetNodeLocations() const;

  // 03 Quadrature
  /**Makes the volumetric/internal quadrature point data for this element.*/
  virtual finite_element::VolumetricQuadraturePointData
  MakeVolumetricQuadraturePointData() const = 0;

  /**Makes the surface quadrature point data for this element, at the specified
   * face.*/
  virtual finite_element::SurfaceQuadraturePointData
  MakeSurfaceQuadraturePointData(size_t face_index) const = 0;

  virtual ~CellMapping() = default;

protected:
  /**This function gets called to compute the cell-volume and
   * face-areas. If simple linear cells are used then the
   * default CellMapping::ComputeCellVolumeAndAreas can be
   * used as this function. Otherwise (i.e. for higher order
   * elements, the child-class should
   * bind a different function to this.*/
  typedef std::function<void(const chi_mesh::MeshContinuum&,
                             const chi_mesh::Cell&,
                             double&,
                             std::vector<double>&)>
    VandAFunction;

  CellMapping(const chi_mesh::MeshContinuum& grid,
              const chi_mesh::Cell& cell,
              size_t num_nodes,
              std::vector<chi_mesh::Vector3> node_locations,
              std::vector<std::vector<int>> face_node_mappings,
              const VandAFunction& volume_area_function);

  /**Static method that all child elements can use as a default.*/
  static void ComputeCellVolumeAndAreas(const chi_mesh::MeshContinuum& grid,
                                        const chi_mesh::Cell& cell,
                                        double& volume,
                                        std::vector<double>& areas);

  const chi_mesh::MeshContinuum& ref_grid_;
  const chi_mesh::Cell& cell_;

  const size_t num_nodes_;
  const std::vector<chi_mesh::Vector3> node_locations_;

  double volume_ = 0.0;
  std::vector<double> areas_;

  /** For each cell face, map from the face node index to the corresponding
   *  cell node index. More specifically, \p face_dof_mappings[f][fi], with
   *  \p fi the face node index of the face identified by face index \p f,
   *  contains the corresponding cell node index. */
  const std::vector<std::vector<int>> face_node_mappings_;
};
} // namespace chi_math

#endif // CHITECH_CELLMAPPING_H
