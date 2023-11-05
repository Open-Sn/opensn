#include "opensn/framework/math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearPolyhedronMapping.h"
#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "opensn/framework/logging/chi_log.h"
#include "opensn/framework/math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

namespace chi_math::cell_mapping
{

PieceWiseLinearPolyhedronMapping::PieceWiseLinearPolyhedronMapping(
  const chi_mesh::Cell& polyh_cell,
  const chi_mesh::MeshContinuum& ref_grid,
  const chi_math::QuadratureTetrahedron& volume_quadrature,
  const chi_math::QuadratureTriangle& surface_quadrature)
  : PieceWiseLinearBaseMapping(ref_grid,
                               polyh_cell,
                               polyh_cell.vertex_ids_.size(), // num_nodes
                               MakeFaceNodeMapping(polyh_cell)),
    volume_quadrature_(volume_quadrature),
    surface_quadrature_(surface_quadrature)
{
  //=========================================== Assign cell centre
  const chi_mesh::Vertex& vcc = polyh_cell.centroid_;
  alphac_ = 1.0 / static_cast<double>(polyh_cell.vertex_ids_.size());

  //=========================================== For each face
  size_t num_faces = polyh_cell.faces_.size();
  face_data_.reserve(num_faces);
  face_betaf_.reserve(num_faces);
  for (size_t f = 0; f < num_faces; f++)
  {
    const chi_mesh::CellFace& face = polyh_cell.faces_[f];
    FEface_data face_f_data;

    face_f_data.normal = face.normal_;

    face_betaf_.push_back(1.0 / static_cast<double>(face.vertex_ids_.size()));

    const chi_mesh::Vertex& vfc = face.centroid_;

    //==================================== For each edge
    const size_t num_edges = face.vertex_ids_.size();
    face_f_data.sides.reserve(num_edges);
    for (size_t e = 0; e < num_edges; ++e)
    {
      FEside_data3d side_data;

      //============================= Assign vertices of tetrahedron
      size_t ep1 = (e < (num_edges - 1)) ? e + 1 : 0;
      uint64_t v0index = face.vertex_ids_[e];
      uint64_t v1index = face.vertex_ids_[ep1];
      side_data.v_index.resize(2, -1);
      side_data.v_index[0] = v0index;
      side_data.v_index[1] = v1index;

      const auto& v0 = ref_grid_.vertices[v0index];
      const auto& v1 = vfc;
      const auto& v2 = ref_grid_.vertices[v1index];
      const auto& v3 = vcc;

      side_data.v0 = v0;

      //============================= Compute vectors
      chi_mesh::Vector3 v01 = v1 - v0;
      chi_mesh::Vector3 v02 = v2 - v0;
      chi_mesh::Vector3 v03 = v3 - v0;

      //============================= Compute determinant of surface jacobian
      // First we compute the rotation matrix which will rotate
      // any vector in natural coordinates to the same reference
      // frame as the current face.
      chi_mesh::Vector3 normal = face.normal_ * -1.0;
      chi_mesh::Vector3 tangent = v02.Cross(normal);
      tangent = tangent / tangent.Norm();
      chi_mesh::Vector3 binorm = v02 / v02.Norm();

      chi_mesh::Matrix3x3 R;
      R.SetColJVec(0, tangent);
      R.SetColJVec(1, binorm);
      R.SetColJVec(2, normal);

      // Now we compute the inverse of this matrix which
      // will allow us to rotate any vector in the same reference
      // frame as the face, to natural coordinates
      chi_mesh::Matrix3x3 Rinv = R.Inverse();

      // Compute v01 and v02 rotated to natural coordinates
      // A test to see if this is done correctly would be to
      // check if fabs(v01N.z) < epsilon and fabs(v02N.z) < epsilon
      chi_mesh::Vector3 v01N = Rinv * v01;
      chi_mesh::Vector3 v02N = Rinv * v02;
      side_data.detJ_surf = v01N.x * v02N.y - v01N.y * v02N.x;

      //============================= Compute Jacobian
      chi_mesh::Matrix3x3 J;
      J.SetColJVec(0, v01);
      J.SetColJVec(1, v02);
      J.SetColJVec(2, v03);

      side_data.J = J;

      //============================= Compute determinant of jacobian
      side_data.detJ = J.Det();

      //============================= Compute inverse Jacobian elements
      chi_mesh::Matrix3x3 JT = J.Transpose();
      chi_mesh::Matrix3x3 Jinv = J.Inverse();
      chi_mesh::Matrix3x3 JTinv = JT.Inverse();

      side_data.Jinv = Jinv;
      side_data.JTinv = JTinv;

      face_f_data.sides.push_back(side_data);
    } // for each edge

    face_data_.push_back(face_f_data);
  } // for each face

  //================================================ Compute Node-Face-Side
  //                                                 mapping
  // This section determines the scope of dof_i on
  // each side (tet) of the cell. If dof_i is on
  // either of the primary tet nodes, it is given
  // index 0 or 1 (-1 otherwise). If the index is
  // -1 the corresponding shapefunction will be
  // ignored when using Ni. The flag "part_of_face"
  // is set when dof_i is part of the same face to
  // which the side belongs and consequently allows
  // the determination of Nf. Nc is always evaluated
  // so no mapping is needed.
  for (int i = 0; i < num_nodes_; i++)
  {
    FEnodeMap newNodeMap;
    for (size_t f = 0; f < face_data_.size(); f++)
    {
      FEnodeFaceMap newFaceMap;
      for (size_t s = 0; s < face_data_[f].sides.size(); s++)
      {
        FEnodeSideMap newSideMap;
        newSideMap.part_of_face = false;
        const uint64_t s0 = face_data_[f].sides[s].v_index[0];
        const uint64_t s1 = face_data_[f].sides[s].v_index[1];
        if (polyh_cell.vertex_ids_[i] == s0)
        {
          newSideMap.index = 0;
          newSideMap.part_of_face = true;
        }
        else if (polyh_cell.vertex_ids_[i] == s1)
        {
          newSideMap.index = 2;
          newSideMap.part_of_face = true;
        }
        else
        {
          newSideMap.index = -1;
          for (size_t v = 0; v < polyh_cell.faces_[f].vertex_ids_.size(); v++)
          {
            if (polyh_cell.vertex_ids_[i] == polyh_cell.faces_[f].vertex_ids_[v])
            {
              newSideMap.part_of_face = true;
              break;
            }
          }
        }
        newFaceMap.side_map.push_back(newSideMap);
      } // for s
      newNodeMap.face_map.push_back(newFaceMap);
    } // for f
    node_side_maps_.push_back(newNodeMap);
  } // for i
}

double
PieceWiseLinearPolyhedronMapping::TetShape(uint32_t index,
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

double
PieceWiseLinearPolyhedronMapping::TetGradShape_x(const uint32_t index)
{
  double value = 0.0;
  if (index == 0) { value = -1.0; }
  if (index == 1) { value = 1.0; }
  if (index == 2) { value = 0.0; }
  if (index == 3) { value = 0.0; }

  return value;
}

double
PieceWiseLinearPolyhedronMapping::TetGradShape_y(const uint32_t index)
{
  double value = 0.0;
  if (index == 0) { value = -1.0; }
  if (index == 1) { value = 0.0; }
  if (index == 2) { value = 1.0; }
  if (index == 3) { value = 0.0; }

  return value;
}

double
PieceWiseLinearPolyhedronMapping::TetGradShape_z(const uint32_t index)
{
  double value = 0.0;
  if (index == 0) { value = -1.0; }
  if (index == 1) { value = 0.0; }
  if (index == 2) { value = 0.0; }
  if (index == 3) { value = 1.0; }

  return value;
}

double
PieceWiseLinearPolyhedronMapping::FaceSideShape(uint32_t face_index,
                                                uint32_t side_index,
                                                uint32_t i,
                                                const chi_mesh::Vector3& qpoint,
                                                bool on_surface /*=false*/) const
{
  double value = 0.0;
  int index = node_side_maps_[i].face_map[face_index].side_map[side_index].index;
  double betaf = face_betaf_[face_index];

  value += TetShape(index, qpoint, on_surface);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    value += betaf * TetShape(1, qpoint, on_surface);
  }
  value += alphac_ * TetShape(3, qpoint, on_surface);

  return value;
}

double
PieceWiseLinearPolyhedronMapping::FaceSideGradShape_x(uint32_t face_index,
                                                      uint32_t side_index,
                                                      uint32_t i) const
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int index = node_side_maps_[i].face_map[face_index].side_map[side_index].index;
  double betaf = face_betaf_[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdx += betaf * TetGradShape_x(1);
  }
  tetdfdx += alphac_ * TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdy += betaf * TetGradShape_y(1);
  }
  tetdfdy += alphac_ * TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdz += betaf * TetGradShape_z(1);
  }
  tetdfdz += alphac_ * TetGradShape_z(3);

  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(0, 0) * tetdfdx;
  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(0, 1) * tetdfdy;
  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(0, 2) * tetdfdz;

  return value;
}

double
PieceWiseLinearPolyhedronMapping::FaceSideGradShape_y(uint32_t face_index,
                                                      uint32_t side_index,
                                                      uint32_t i) const
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int index = node_side_maps_[i].face_map[face_index].side_map[side_index].index;
  double betaf = face_betaf_[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdx += betaf * TetGradShape_x(1);
  }
  tetdfdx += alphac_ * TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdy += betaf * TetGradShape_y(1);
  }
  tetdfdy += alphac_ * TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdz += betaf * TetGradShape_z(1);
  }
  tetdfdz += alphac_ * TetGradShape_z(3);

  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(1, 0) * tetdfdx;
  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(1, 1) * tetdfdy;
  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(1, 2) * tetdfdz;

  return value;
}

double
PieceWiseLinearPolyhedronMapping::FaceSideGradShape_z(uint32_t face_index,
                                                      uint32_t side_index,
                                                      uint32_t i) const
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int index = node_side_maps_[i].face_map[face_index].side_map[side_index].index;
  double betaf = face_betaf_[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdx += betaf * TetGradShape_x(1);
  }
  tetdfdx += alphac_ * TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdy += betaf * TetGradShape_y(1);
  }
  tetdfdy += alphac_ * TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps_[i].face_map[face_index].side_map[side_index].part_of_face)
  {
    tetdfdz += betaf * TetGradShape_z(1);
  }
  tetdfdz += alphac_ * TetGradShape_z(3);

  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(2, 0) * tetdfdx;
  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(2, 1) * tetdfdy;
  value += face_data_[face_index].sides[side_index].JTinv.GetIJ(2, 2) * tetdfdz;

  return value;
}

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

void
PieceWiseLinearPolyhedronMapping::GradShapeValues(
  const chi_mesh::Vector3& xyz, std::vector<chi_mesh::Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    gradshape_values.emplace_back(GradShapeValue(i, xyz));
}

finite_element::VolumetricQuadraturePointData
PieceWiseLinearPolyhedronMapping::MakeVolumetricQuadraturePointData() const
{
  //=================================== Determine number of internal qpoints
  size_t num_tets = 0;
  for (auto& face : face_data_)
    num_tets += face.sides.size();

  size_t num_vol_qpoints = volume_quadrature_.qpoints_.size();
  size_t ttl_num_vol_qpoints = num_tets * num_vol_qpoints;

  //=================================== Declare necessary vars
  std::vector<unsigned int> V_quadrature_point_indices;
  VecVec3 V_qpoints_xyz;
  std::vector<VecDbl> V_shape_value;
  std::vector<VecVec3> V_shape_grad;
  VecDbl V_JxW;
  size_t V_num_nodes;

  //=================================== Init volumetric quadrature
  V_quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp = 0; qp < ttl_num_vol_qpoints; ++qp)
    V_quadrature_point_indices.push_back(qp);

  V_shape_value.reserve(num_nodes_);
  V_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; i++)
  {
    VecDbl node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (size_t f = 0; f < face_data_.size(); f++)
    {
      for (size_t s = 0; s < face_data_[f].sides.size(); s++)
      {
        for (const auto& qpoint : volume_quadrature_.qpoints_)
        {
          node_shape_value.push_back(FaceSideShape(f, s, i, qpoint));
          node_shape_grad.emplace_back(FaceSideGradShape_x(f, s, i),  // x
                                       FaceSideGradShape_y(f, s, i),  // y
                                       FaceSideGradShape_z(f, s, i)); // z
        }                                                             // for qp
      }                                                               // for side
    }                                                                 // for face

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  } // for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (const auto& face : face_data_)
  {
    for (const auto& side : face.sides)
    {
      for (size_t qp = 0; qp < num_vol_qpoints; ++qp)
      {
        const auto w = volume_quadrature_.weights_[qp];
        V_JxW.push_back(side.detJ * w);

        const auto& qp_xyz_tilde = volume_quadrature_.qpoints_[qp];
        V_qpoints_xyz.push_back(side.v0 + side.J * qp_xyz_tilde);
      } // for qp
    }   // for side
  }     // for face

  V_num_nodes = num_nodes_;

  return finite_element::VolumetricQuadraturePointData(V_quadrature_point_indices,
                                                       V_qpoints_xyz,
                                                       V_shape_value,
                                                       V_shape_grad,
                                                       V_JxW,
                                                       face_node_mappings_,
                                                       V_num_nodes);
}

finite_element::SurfaceQuadraturePointData
PieceWiseLinearPolyhedronMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  const bool ON_SURFACE = true;

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = surface_quadrature_.qpoints_.size();

  unsigned int f = face_index;
  //=================================== Declare necessary vars
  std::vector<unsigned int> F_quadrature_point_indices;
  VecVec3 F_qpoints_xyz;
  std::vector<VecDbl> F_shape_value;
  std::vector<VecVec3> F_shape_grad;
  VecDbl F_JxW;
  VecVec3 F_normals;
  size_t F_num_nodes;

  size_t num_tris = face_data_[f].sides.size();
  size_t ttl_num_face_qpoints = num_tris * num_srf_qpoints;

  F_quadrature_point_indices.reserve(ttl_num_face_qpoints);
  for (unsigned int qp = 0; qp < ttl_num_face_qpoints; ++qp)
    F_quadrature_point_indices.push_back(qp);

  F_normals.reserve(ttl_num_face_qpoints);
  for (size_t qp = 0; qp < ttl_num_face_qpoints; ++qp)
    F_normals.push_back(face_data_[f].normal);

  F_shape_value.reserve(num_nodes_);
  F_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; i++)
  {
    VecDbl node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_face_qpoints);
    node_shape_grad.reserve(ttl_num_face_qpoints);

    for (size_t s = 0; s < face_data_[f].sides.size(); s++)
    {
      for (const auto& qpoint : surface_quadrature_.qpoints_)
      {
        node_shape_value.push_back(FaceSideShape(f, s, i, qpoint, ON_SURFACE));
        node_shape_grad.emplace_back(FaceSideGradShape_x(f, s, i),  // x
                                     FaceSideGradShape_y(f, s, i),  // y
                                     FaceSideGradShape_z(f, s, i)); // z
      }                                                             // for qp
    }                                                               // for s
    F_shape_value.push_back(node_shape_value);
    F_shape_grad.push_back(node_shape_grad);
  } // for i

  F_JxW.reserve(ttl_num_face_qpoints);
  F_qpoints_xyz.reserve(ttl_num_face_qpoints);
  for (const auto& side : face_data_[f].sides)
    for (size_t qp = 0; qp < num_srf_qpoints; ++qp)
    {
      const auto w = surface_quadrature_.weights_[qp];
      F_JxW.push_back(side.detJ_surf * w);

      const auto& qp_xyz_tilde = surface_quadrature_.qpoints_[qp];
      F_qpoints_xyz.push_back(side.v0 + side.J * qp_xyz_tilde);
    }

  F_num_nodes = face_data_[f].sides.size();

  return finite_element::SurfaceQuadraturePointData(F_quadrature_point_indices,
                                                    F_qpoints_xyz,
                                                    F_shape_value,
                                                    F_shape_grad,
                                                    F_JxW,
                                                    F_normals,
                                                    face_node_mappings_,
                                                    F_num_nodes);
}

} // namespace chi_math::cell_mapping
