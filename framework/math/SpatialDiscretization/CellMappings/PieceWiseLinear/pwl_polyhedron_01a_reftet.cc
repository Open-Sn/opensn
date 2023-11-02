#include "PieceWiseLinearPolyhedronMapping.h"

namespace chi_math::cell_mapping
{

double PieceWiseLinearPolyhedronMapping::TetShape(uint32_t index,
                                         const chi_mesh::Vector3& qpoint,
                                         bool on_surface /*=false*/)
{
  double value = 0.0;

  if (index == 0) { value = 1.0 - qpoint.x - qpoint.y - qpoint.z; }
  if (index == 1) { value = qpoint.x; }
  if (index == 2) { value = qpoint.y; }
  if (index == 3) { value = qpoint.z; }

  return value;
}

double PieceWiseLinearPolyhedronMapping::TetGradShape_x(const uint32_t index)
{
  double value = 0.0;
  if (index == 0) { value = -1.0; }
  if (index == 1) { value = 1.0; }
  if (index == 2) { value = 0.0; }
  if (index == 3) { value = 0.0; }

  return value;
}

double PieceWiseLinearPolyhedronMapping::TetGradShape_y(const uint32_t index)
{
  double value = 0.0;
  if (index == 0) { value = -1.0; }
  if (index == 1) { value = 0.0; }
  if (index == 2) { value = 1.0; }
  if (index == 3) { value = 0.0; }

  return value;
}

double PieceWiseLinearPolyhedronMapping::TetGradShape_z(const uint32_t index)
{
  double value = 0.0;
  if (index == 0) { value = -1.0; }
  if (index == 1) { value = 0.0; }
  if (index == 2) { value = 0.0; }
  if (index == 3) { value = 1.0; }

  return value;
}

} // namespace chi_math::cell_mapping
