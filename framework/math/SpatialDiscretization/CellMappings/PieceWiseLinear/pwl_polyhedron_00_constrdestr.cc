#include "PieceWiseLinearPolyhedronMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"

namespace chi_math::cell_mapping
{

// ###################################################################
/**Constructor for the Piecewise Linear Polyhedron cell finite elment
 * view.
 *
 * */
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
            if (polyh_cell.vertex_ids_[i] ==
                polyh_cell.faces_[f].vertex_ids_[v])
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

} // namespace chi_math::cell_mapping
