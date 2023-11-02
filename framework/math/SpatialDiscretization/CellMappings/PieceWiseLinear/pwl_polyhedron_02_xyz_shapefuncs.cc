#include "PieceWiseLinearPolyhedronMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_math::cell_mapping
{

// ###################################################################
/**Returns the evaluation of shape function i at the supplied point.*/
double
PieceWiseLinearPolyhedronMapping::ShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  for (size_t f = 0; f < face_data_.size(); f++)
  {
    for (size_t s = 0; s < face_data_[f].sides.size(); s++)
    {
      // Map xyz to xi_eta_zeta
      const auto& p0 = ref_grid_.vertices[face_data_[f].sides[s].v_index[0]];
      chi_mesh::Vector3 xyz_ref = xyz - p0;

      chi_mesh::Vector3 xi_eta_zeta = face_data_[f].sides[s].Jinv * xyz_ref;

      double xi = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta = xi_eta_zeta.z;

      // Determine if inside tet
      if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and (zeta >= -1.0e-12) and
          ((xi + eta + zeta) <= (1.0 + 1.0e-12)))
      {
        double Ni = 0.0;
        double Nf = 0.0;
        double Nc = alphac_ * zeta;

        if (node_side_maps_[i].face_map[f].side_map[s].part_of_face)
        {
          if (node_side_maps_[i].face_map[f].side_map[s].index == 0) { Ni = 1 - xi - eta - zeta; }
          if (node_side_maps_[i].face_map[f].side_map[s].index == 2) { Ni = eta; }

          Nf = face_betaf_[f] * xi;
        }

        return Ni + Nf + Nc;
      }
    }
  }
  return 0.0;
}

// ###################################################################
/**Populates shape_values with the value of each shape function's
 * value evaluate at the supplied point.*/
void
PieceWiseLinearPolyhedronMapping::ShapeValues(const chi_mesh::Vector3& xyz,
                                              std::vector<double>& shape_values) const
{
  shape_values.resize(num_nodes_, 0.0);
  for (size_t f = 0; f < face_data_.size(); f++)
  {
    for (size_t s = 0; s < face_data_[f].sides.size(); s++)
    {
      auto& side_fe_info = face_data_[f].sides[s];
      // Map xyz to xi_eta_zeta
      const auto& p0 = ref_grid_.vertices[side_fe_info.v_index[0]];
      chi_mesh::Vector3 xi_eta_zeta = side_fe_info.Jinv * (xyz - p0);

      double xi = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta = xi_eta_zeta.z;

      // Determine if inside tet
      if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and (zeta >= -1.0e-12) and
          ((xi + eta + zeta) <= (1.0 + 1.0e-12)))
      {
        for (int i = 0; i < num_nodes_; i++)
        {
          auto side_map = node_side_maps_[i].face_map[f].side_map[s];

          double Ni = 0.0;
          double Nf = 0.0;
          double Nc = alphac_ * zeta;

          if (side_map.part_of_face)
          {
            if (side_map.index == 0) Ni = 1 - xi - eta - zeta;
            else if (side_map.index == 2)
              Ni = eta;

            Nf = face_betaf_[f] * xi;
          }

          shape_values[i] = Ni + Nf + Nc;
        } // for dof
        return;
      } // if in tet
    }   // for side
  }     // for face
}

// ###################################################################
/**Returns the evaluation of grad-shape function i at the supplied point.*/
chi_mesh::Vector3
PieceWiseLinearPolyhedronMapping::GradShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  chi_mesh::Vector3 grad, gradr;
  for (size_t f = 0; f < face_data_.size(); f++)
  {
    for (size_t s = 0; s < face_data_[f].sides.size(); s++)
    {
      // Map xyz to xi_eta_zeta
      const auto& p0 = ref_grid_.vertices[face_data_[f].sides[s].v_index[0]];
      chi_mesh::Vector3 xyz_ref = xyz - p0;

      chi_mesh::Vector3 xi_eta_zeta = face_data_[f].sides[s].Jinv * xyz_ref;

      double xi = xi_eta_zeta.x;
      double eta = xi_eta_zeta.y;
      double zeta = xi_eta_zeta.z;

      // Determine if inside tet
      if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and (zeta >= -1.0e-12) and
          ((xi + eta + zeta) <= (1.0 + 1.0e-12)))
      {
        chi_mesh::Vector3 grad_i;
        chi_mesh::Vector3 grad_f;
        chi_mesh::Vector3 grad_c;

        if (node_side_maps_[i].face_map[f].side_map[s].part_of_face)
        {
          if (node_side_maps_[i].face_map[f].side_map[s].index == 0)
          {
            grad_i.x = -1.0;
            grad_i.y = -1.0;
            grad_i.z = -1.0;
          }
          if (node_side_maps_[i].face_map[f].side_map[s].index == 2)
          {
            grad_i.x = 0.0;
            grad_i.y = 1.0;
            grad_i.z = 0.0;
          }

          grad_f.x = face_betaf_[f] * 1.0;
          grad_f.y = face_betaf_[f] * 0.0;
          grad_f.z = face_betaf_[f] * 0.0;
        }

        grad_c.x = alphac_ * 0.0;
        grad_c.y = alphac_ * 0.0;
        grad_c.z = alphac_ * 1.0;

        grad = (grad_i + grad_f + grad_c);
        grad = face_data_[f].sides[s].JTinv * grad;

        return grad;
      }
    }
  }
  return gradr;
}

// ###################################################################
/**Populates gradshape_values with the value of each shape function's
 * gradient evaluated at the supplied point.*/
void
PieceWiseLinearPolyhedronMapping::GradShapeValues(
  const chi_mesh::Vector3& xyz, std::vector<chi_mesh::Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    gradshape_values.emplace_back(GradShapeValue(i, xyz));
}

} // namespace chi_math::cell_mapping
