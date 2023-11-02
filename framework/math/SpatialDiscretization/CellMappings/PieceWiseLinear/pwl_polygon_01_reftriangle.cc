#include "PieceWiseLinearPolygonMapping.h"

namespace chi_math::cell_mapping
{

double
PieceWiseLinearPolygonMapping::TriShape(uint32_t index,
                                        const chi_mesh::Vector3& qpoint,
                                        bool on_surface /*false*/)
{
  double xi;
  double eta;
  if (!on_surface)
  {
    xi = qpoint.x;
    eta = qpoint.y;
  }
  else
  {
    xi = qpoint.x;
    eta = 0.0;
  }

  double value = 0.0;
  if (index == 0) value = 1.0 - xi - eta;
  else if (index == 1)
    value = xi;
  else if (index == 2)
    value = eta;

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Varphi_x
/**Precomputation of the shape function at a quadrature point.*/
double
PieceWiseLinearPolygonMapping::SideShape(uint32_t side,
                                         uint32_t i,
                                         const chi_mesh::Vector3& qpoint,
                                         bool on_surface /*=false*/) const
{
  int index = node_to_side_map_[i][side];
  double value = 0.0;
  if (index == 0 or index == 1) value = TriShape(index, qpoint, on_surface);

  value += beta_ * TriShape(2, qpoint, on_surface);

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_x
/**Precomputation of the partial derivative along x of the
 * shape function at a quadrature point.*/
double
PieceWiseLinearPolygonMapping::SideGradShape_x(uint32_t side, uint32_t i) const
{
  int index = node_to_side_map_[i][side];
  double value = 0;
  if (index == 0)
  {

    value = sides_[side].JTinv.GetIJ(0, 0) * -1.0 + sides_[side].JTinv.GetIJ(0, 1) * -1.0;
  }
  if (index == 1)
  {

    value = sides_[side].JTinv.GetIJ(0, 0) * 1.0 + sides_[side].JTinv.GetIJ(0, 1) * 0.0;
  }

  value += beta_ * (sides_[side].JTinv.GetIJ(0, 0) * 0.0 + sides_[side].JTinv.GetIJ(0, 1) * 1.0);

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_y
/**Precomputation of the partial derivative along y of the
 * shape function at a quadrature point.*/
double
PieceWiseLinearPolygonMapping::SideGradShape_y(uint32_t side, uint32_t i) const
{
  int index = node_to_side_map_[i][side];
  double value = 0;
  if (index == 0)
  {

    value = sides_[side].JTinv.GetIJ(1, 0) * -1.0 + sides_[side].JTinv.GetIJ(1, 1) * -1.0;
  }
  if (index == 1)
  {

    value = sides_[side].JTinv.GetIJ(1, 0) * 1.0 + sides_[side].JTinv.GetIJ(1, 1) * 0.0;
  }

  value += beta_ * (sides_[side].JTinv.GetIJ(1, 0) * 0.0 + sides_[side].JTinv.GetIJ(1, 1) * 1.0);

  return value;
}

} // namespace chi_math::cell_mapping
