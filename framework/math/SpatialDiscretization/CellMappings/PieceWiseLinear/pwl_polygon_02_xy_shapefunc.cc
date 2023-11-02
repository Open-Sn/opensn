#include "PieceWiseLinearPolygonMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_math::cell_mapping
{

//###################################################################
/**Returns the value of the shape function given cartesian
 * coordinates.*/
double
PieceWiseLinearPolygonMapping::ShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  for (int s = 0; s < num_of_subtris_; s++)
  {
    const auto& p0 = ref_grid_.vertices[sides_[s].v_index[0]];
    chi_mesh::Vector3 xyz_ref = xyz - p0;

    chi_mesh::Vector3 xi_eta_zeta = sides_[s].Jinv * xyz_ref;

    double xi = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    // Determine if inside tet
    if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and ((xi + eta) <= (1.0 + 1.0e-12)))
    {
      int index = node_to_side_map_[i][s];
      double value = 0.0;

      if (index == 0) { value = 1.0 - xi - eta; }
      if (index == 1) { value = xi; }

      value += beta_ * eta;

      return value;
    }
  }

  return 0.0;
}

//###################################################################
/**Populates shape_values with the value of each shape function's
 * value evaluate at the supplied point.*/
void
PieceWiseLinearPolygonMapping::ShapeValues(const chi_mesh::Vector3& xyz,
                                           std::vector<double>& shape_values) const
{
  shape_values.resize(num_nodes_, 0.0);
  for (int s = 0; s < num_of_subtris_; s++)
  {
    const auto& p0 = ref_grid_.vertices[sides_[s].v_index[0]];
    chi_mesh::Vector3 xi_eta_zeta = sides_[s].Jinv * (xyz - p0);

    double xi = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    // Determine if inside tet
    if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and ((xi + eta) <= (1.0 + 1.0e-12)))
    {
      for (int i = 0; i < num_nodes_; i++)
      {
        int index = node_to_side_map_[i][s];
        double value = 0.0;

        if (index == 0) { value = 1.0 - xi - eta; }
        if (index == 1) { value = xi; }

        value += beta_ * eta;

        shape_values[i] = value;
      }
      return;
    } // if in triangle
  }   // for side
}

//###################################################################
/**Returns the evaluation of grad-shape function i at the supplied point.*/
chi_mesh::Vector3
PieceWiseLinearPolygonMapping::GradShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  chi_mesh::Vector3 grad_r;
  chi_mesh::Vector3 grad;

  for (int e = 0; e < num_of_subtris_; e++)
  {
    const auto& p0 = ref_grid_.vertices[sides_[e].v_index[0]];
    chi_mesh::Vector3 xyz_ref = xyz - p0;

    chi_mesh::Vector3 xi_eta_zeta = sides_[e].Jinv * xyz_ref;

    double xi = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and ((xi + eta) <= (1.0 + 1.0e-12)))
    {
      int index = node_to_side_map_[i][e];

      if (index == 0)
      {
        grad_r.x += -1.0;
        grad_r.y += -1.0;
      }
      if (index == 1)
      {
        grad_r.x += 1.0;
        grad_r.y += 0.0;
      }

      grad_r.y += beta_ * 1.0;

      grad = sides_[e].JTinv * grad_r;

      return grad;
    }
  }

  return grad;
}

//###################################################################
/**Populates gradshape_values with the value of each shape function's
 * gradient evaluated at the supplied point.*/
void
PieceWiseLinearPolygonMapping::GradShapeValues(
  const chi_mesh::Vector3& xyz, std::vector<chi_mesh::Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    gradshape_values.emplace_back(GradShapeValue(i, xyz));
}

} // namespace chi_math::cell_mapping
