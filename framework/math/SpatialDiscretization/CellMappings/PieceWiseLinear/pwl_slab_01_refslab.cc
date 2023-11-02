#include "PieceWiseLinearSlabMapping.h"

namespace chi_math::cell_mapping
{

double PieceWiseLinearSlabMapping::SlabShape(uint32_t index,
                                    const chi_mesh::Vector3& qpoint,
                                    bool on_surface,
                                    uint32_t edge) const
{
  double xi = 0.0;
  if (!on_surface) xi = qpoint.x;
  else
    xi = static_cast<double>(edge);

  double value = 0.0;
  if (index == 0) value = 1.0 - xi;
  else if (index == 1)
    value = xi;

  return value;
}

double PieceWiseLinearSlabMapping::SlabGradShape(uint32_t index) const
{
  double value = 0.0;

  if (index == 0) value = -1.0 / h_;
  else
    value = 1.0 / h_;

  return value;
}

} // namespace chi_math::cell_mapping
