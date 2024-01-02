#include "framework/math/spatial_discretization/cell_mappings/lagrange_base_mapping.h"

#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/math/spatial_discretization/finite_element/quadrature_point_data.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

LagrangeBaseMapping::LagrangeBaseMapping(const MeshContinuum& grid,
                                         const Cell& cell,
                                         size_t num_nodes,
                                         std::vector<std::vector<int>> face_node_mappings,
                                         const Quadrature& volume_quadrature,
                                         const Quadrature& surface_quadrature)
  : CellMapping(grid,
                cell,
                num_nodes,
                GetVertexLocations(grid, cell),
                std::move(face_node_mappings),
                &CellMapping::ComputeCellVolumeAndAreas),
    volume_quadrature_(volume_quadrature),
    surface_quadrature_(surface_quadrature)
{
}

std::vector<std::vector<int>>
LagrangeBaseMapping::MakeFaceNodeMapping(const Cell& cell)
{
  const size_t num_faces = cell.faces_.size();
  std::vector<std::vector<int>> mappings;
  mappings.reserve(num_faces);
  for (auto& face : cell.faces_)
  {
    std::vector<int> face_dof_mapping;
    face_dof_mapping.reserve(face.vertex_ids_.size());
    for (uint64_t fvid : face.vertex_ids_)
    {
      int mapping = -1;
      for (size_t ci = 0; ci < cell.vertex_ids_.size(); ci++)
      {
        if (fvid == cell.vertex_ids_[ci])
        {
          mapping = static_cast<int>(ci);
          break;
        }
      } // for cell i
      if (mapping < 0)
      {
        log.LogAllError() << "Unknown face mapping encountered. "
                             "pwl_polyhedron.h";
        Exit(EXIT_FAILURE);
      }
      face_dof_mapping.push_back(mapping);
    } // for face i

    mappings.push_back(face_dof_mapping);
  }
  return mappings;
}

std::vector<Vector3>
LagrangeBaseMapping::GetVertexLocations(const MeshContinuum& grid, const Cell& cell)
{
  std::vector<Vector3> verts;
  verts.reserve(cell.vertex_ids_.size());

  for (const auto vid : cell.vertex_ids_)
    verts.push_back(grid.vertices[vid]);

  return verts;
}

Vector3
LagrangeBaseMapping::MapWorldXYZToNaturalXYZ(const Vec3& world_xyz) const
{
  WorldXYZToNaturalMappingHelper helper(*this, world_xyz);

  auto nat_xyz = NewtonIteration(helper, {0.0, 0.0, 0.0}, 20, 1.0e-8);

  return {nat_xyz[0], nat_xyz[1], nat_xyz[2]};
}

double
LagrangeBaseMapping::ShapeValue(int i, const Vector3& xyz) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);
  return RefShape(i, natural_xyz);
}

void
LagrangeBaseMapping::ShapeValues(const Vector3& xyz, std::vector<double>& shape_values) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);

  shape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    shape_values.push_back(RefShape(i, natural_xyz));
}

Vector3
LagrangeBaseMapping::GradShapeValue(int i, const Vector3& xyz) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);
  return RefGradShape(i, natural_xyz);
}

void
LagrangeBaseMapping::GradShapeValues(const Vector3& xyz,
                                     std::vector<Vector3>& gradshape_values) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);

  gradshape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    gradshape_values.push_back(RefGradShape(i, natural_xyz));
}

VolumetricQuadraturePointData
LagrangeBaseMapping::MakeVolumetricQuadraturePointData() const
{
  typedef std::vector<Vec3> VecVec3;
  const size_t num_qpoints = volume_quadrature_.qpoints_.size();

  // Declare necessary vars
  std::vector<unsigned int> quadrature_point_indices;
  VecVec3 qpoints_xyz;
  std::vector<VecDbl> shape_value;
  std::vector<VecVec3> shape_grad;
  VecDbl JxW;

  quadrature_point_indices.reserve(num_qpoints);
  for (unsigned int qp = 0; qp < num_qpoints; ++qp)
    quadrature_point_indices.push_back(qp);

  JxW.assign(num_qpoints, 0.0);
  shape_value.assign(num_nodes_, VecDbl(num_qpoints, 0.0));
  shape_grad.assign(num_nodes_, VecVec3(num_qpoints));
  qpoints_xyz.assign(num_qpoints, Vec3{0.0, 0.0, 0.0});

  for (uint32_t qp : quadrature_point_indices)
  {
    const Vec3& qpoint = volume_quadrature_.qpoints_[qp];
    const auto J = RefJacobian(qpoint);
    const auto JT = Transpose(J);
    const auto JTinv = Inverse(JT);
    const double detJ = Determinant(J);

    JxW[qp] = detJ * volume_quadrature_.weights_[qp];

    for (size_t i = 0; i < num_nodes_; ++i)
    {
      const double ref_shape_i = RefShape(i, qpoint);
      shape_value[i][qp] = ref_shape_i;
      const Vec3 ref_shape_grad = RefGradShape(i, qpoint);
      const VecDbl b = {ref_shape_grad.x, ref_shape_grad.y, ref_shape_grad.z};
      const VecDbl JTInv_b = MatMul(JTinv, b);

      shape_grad[i][qp] = {JTInv_b[0], JTInv_b[1], JTInv_b[2]};

      const auto& node_xyz = node_locations_[i];

      qpoints_xyz[qp] += ref_shape_i * node_xyz;
    } // for i
  }   // for qp

  return VolumetricQuadraturePointData(quadrature_point_indices,
                                       qpoints_xyz,
                                       shape_value,
                                       shape_grad,
                                       JxW,
                                       face_node_mappings_,
                                       num_nodes_);
}

const Quadrature&
LagrangeBaseMapping::GetSurfaceQuadrature(size_t face_index) const
{
  return surface_quadrature_;
}

SurfaceQuadraturePointData
LagrangeBaseMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  const auto& surface_quadrature = GetSurfaceQuadrature(face_index);

  typedef std::vector<Vec3> VecVec3;
  const size_t num_qpoints = surface_quadrature.qpoints_.size();

  // Declare necessary vars
  std::vector<unsigned int> quadrature_point_indices;
  VecVec3 qpoints_xyz;
  std::vector<VecDbl> shape_value;
  std::vector<VecVec3> shape_grad;
  VecDbl JxW;
  VecVec3 normals;

  quadrature_point_indices.reserve(num_qpoints);
  for (unsigned int qp = 0; qp < num_qpoints; ++qp)
    quadrature_point_indices.push_back(qp);

  JxW.assign(num_qpoints, 0.0);
  normals.assign(num_qpoints, Vec3{0.0, 0.0, 0.0});
  shape_value.assign(num_nodes_, VecDbl(num_qpoints, 0.0));
  shape_grad.assign(num_nodes_, VecVec3(num_qpoints));
  qpoints_xyz.assign(num_qpoints, Vec3{0.0, 0.0, 0.0});

  const size_t f = face_index;
  for (uint32_t qp : quadrature_point_indices)
  {
    const Vec3& qpoint_face = surface_quadrature.qpoints_[qp];
    const auto qpoint = FaceToElementQPointConversion(f, qpoint_face);

    const auto J = RefJacobian(qpoint);
    const auto JT = Transpose(J);
    const auto JTinv = Inverse(JT);
    const auto [detJ, qp_normal] = RefFaceJacobianDeterminantAndNormal(f, qpoint_face);

    normals[qp] = qp_normal;

    JxW[qp] = detJ * surface_quadrature.weights_[qp];

    for (size_t i = 0; i < num_nodes_; ++i)
    {
      const double ref_shape_i = RefShape(i, qpoint);
      shape_value[i][qp] = ref_shape_i;
      const Vec3 ref_shape_grad = RefGradShape(i, qpoint);
      const VecDbl b = {ref_shape_grad.x, ref_shape_grad.y, ref_shape_grad.z};
      const VecDbl JTInv_b = MatMul(JTinv, b);

      shape_grad[i][qp] = {JTInv_b[0], JTInv_b[1], JTInv_b[2]};

      const auto& node_xyz = node_locations_[i];

      qpoints_xyz[qp] += ref_shape_i * node_xyz;
    } // for i
  }   // for qp

  return SurfaceQuadraturePointData(quadrature_point_indices,
                                    qpoints_xyz,
                                    shape_value,
                                    shape_grad,
                                    JxW,
                                    normals,
                                    face_node_mappings_,
                                    num_nodes_);
}

std::pair<double, LagrangeBaseMapping::Vec3>
LagrangeBaseMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                                         const Vec3& qpoint_face) const
{
  ChiLogicalError("Method not implemented");
}

WorldXYZToNaturalMappingHelper::WorldXYZToNaturalMappingHelper(
  const LagrangeBaseMapping& cell_mapping, const Vec3& world_x)
  : cell_mapping_(cell_mapping), world_x_(world_x)
{
  const auto& cell = cell_mapping.ReferenceCell();
  const auto cell_type = cell.Type();
  if (cell_type == CellType::SLAB) dimension_ = 1;
  else if (cell_type == CellType::POLYGON)
    dimension_ = 2;
  else if (cell_type == CellType::POLYHEDRON)
    dimension_ = 3;
  else
    ChiLogicalError("Unsupported cell type.");
}

VecDbl
WorldXYZToNaturalMappingHelper::F(const VecDbl& x) const
{
  const size_t num_nodes = cell_mapping_.NumNodes();
  const auto& node_x = cell_mapping_.GetNodeLocations();

  Vec3 result_vec3 = world_x_;
  if (dimension_ == 1)
    for (uint32_t i = 0; i < num_nodes; ++i)
    {
      const double shape_i_val = cell_mapping_.RefShape(i, {x[0], x[1], x[2]});
      result_vec3 -= Vec3(0.0, 0.0, shape_i_val * node_x[i].z);
    }
  else if (dimension_ == 2)
    for (uint32_t i = 0; i < num_nodes; ++i)
    {
      const double shape_i_val = cell_mapping_.RefShape(i, {x[0], x[1], x[2]});
      result_vec3 -= Vec3(shape_i_val * node_x[i].x, shape_i_val * node_x[i].y, 0.0);
    }
  else if (dimension_ == 3)
    for (uint32_t i = 0; i < num_nodes; ++i)
      result_vec3 -= cell_mapping_.RefShape(i, {x[0], x[1], x[2]}) * node_x[i];

  return {result_vec3.x, result_vec3.y, result_vec3.z};
}

MatDbl
WorldXYZToNaturalMappingHelper::J(const VecDbl& x) const
{
  return cell_mapping_.RefJacobian({x[0], x[1], x[2]});
}

} // namespace opensn
