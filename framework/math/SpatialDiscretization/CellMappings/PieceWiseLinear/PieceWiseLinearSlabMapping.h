#ifndef PWL_SLAB_VALUES_H
#define PWL_SLAB_VALUES_H

#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"
#include "math/Quadratures/quadrature_line.h"
#include "mesh/Cell/cell.h"
#include <array>

// ###################################################################
namespace chi_math::cell_mapping
{

/**Object for handling slab shaped piecewise linear shape functions.
* \ingroup doc_CellMappings*/
class PieceWiseLinearSlabMapping : public PieceWiseLinearBaseMapping
{
public:
  /**Constructor for a slab view.*/
  PieceWiseLinearSlabMapping(const chi_mesh::Cell& slab_cell,
                      const chi_mesh::MeshContinuum& ref_grid,
                      const QuadratureLine& volume_quadrature);

  finite_element::VolumetricQuadraturePointData
  MakeVolumetricQuadraturePointData() const override;

  finite_element::SurfaceQuadraturePointData
  MakeSurfaceQuadraturePointData(size_t face_index) const override;

  // ################################################## Define standard
  //                                                    slab linear shape
  //                                                    functions
  double SlabShape(uint32_t index,
                   const chi_mesh::Vector3& qpoint,
                   bool on_surface = false,
                   uint32_t edge = 0) const;

  double SlabGradShape(uint32_t index) const;

  // ############################################### Actual shape functions
  //                                                 as function of cartesian
  //                                                 coordinates
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) const override;

  chi_mesh::Vector3 GradShapeValue(int i,
                                   const chi_mesh::Vector3& xyz) const override;

  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) const override;

  void GradShapeValues(
    const chi_mesh::Vector3& xyz,
    std::vector<chi_mesh::Vector3>& gradshape_values) const override;

private:
  chi_mesh::Vector3 v0_;
  uint64_t v0i_;
  uint64_t v1i_;
  std::array<chi_mesh::Normal, 2> normals_;
  const QuadratureLine& volume_quadrature_;
  double h_;
};
} // namespace chi_math::cell_mapping
#endif
