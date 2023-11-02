#include "sldfe_sq.h"

#include <algorithm>

//###################################################################
/**Generates the standard points on the reference face.*/
void
chi_math::SimplifiedLDFESQ::Quadrature::GenerateReferenceFaceVertices(
  const chi_mesh::Matrix3x3& rotation_matrix, const chi_mesh::Vector3& translation, int level)
{
  typedef std::vector<chi_mesh::Vertex> VertList;

  int Ns = (level + 1); // Number of subdivisions
  int Np = Ns + 1;      // Number of diagonal points

  chi_math::QuadratureGaussLegendre legendre(QuadratureOrder::THIRTYSECOND);

  //============================================= Generate xy_tilde values
  std::vector<VertList> vertices_xy_tilde_ij;
  vertices_xy_tilde_ij.resize(Np, VertList(Np));
  for (int i = 0; i < Np; ++i)
    for (int j = 0; j < Np; ++j)
      vertices_xy_tilde_ij[i][j] =
        chi_mesh::Vertex(diagonal_vertices_[i].x, diagonal_vertices_[j].y, 0.0);

  //============================================= Generate SQs
  for (int i = 0; i < Ns; ++i)
  {
    for (int j = 0; j < Ns; ++j)
    {
      SphericalQuadrilateral sq;

      sq.rotation_matrix = rotation_matrix;
      sq.translation_vector = translation;

      //==================================== Set xy-tilde vertices
      sq.vertices_xy_tilde[0] = vertices_xy_tilde_ij[i][j];
      sq.vertices_xy_tilde[1] = vertices_xy_tilde_ij[i + 1][j];
      sq.vertices_xy_tilde[2] = vertices_xy_tilde_ij[i + 1][j + 1];
      sq.vertices_xy_tilde[3] = vertices_xy_tilde_ij[i][j + 1];
      auto& vxy = sq.vertices_xy_tilde;

      //==================================== Set xyz_prime vertices
      for (int v = 0; v < 4; ++v)
        sq.vertices_xyz_prime[v] = rotation_matrix * vxy[v] + translation;

      //==================================== Set xyz vertices
      for (int v = 0; v < 4; ++v)
        sq.vertices_xyz[v] = sq.vertices_xyz_prime[v].Normalized();

      //==================================== Compute SQ xyz-centroid
      for (auto& vertex : sq.vertices_xyz)
        sq.centroid_xyz += vertex;
      sq.centroid_xyz /= 4;
      sq.centroid_xyz.Normalize();

      auto v0 = sq.centroid_xyz.Normalized();
      auto v1 = sq.vertices_xyz[0];
      auto v2 = sq.vertices_xyz[1];

      //==================================== Correction orientation
      if ((v1 - v0).Cross(v2 - v0).Dot(v0) < 0.0)
      {
        std::reverse(sq.vertices_xy_tilde.begin(), sq.vertices_xy_tilde.end());
        std::reverse(sq.vertices_xyz_prime.begin(), sq.vertices_xyz_prime.end());
        std::reverse(sq.vertices_xyz.begin(), sq.vertices_xyz.end());
      }

      //==================================== Compute area
      sq.area = ComputeSphericalQuadrilateralArea(sq.vertices_xyz);

      //==================================== Set octant modifier
      sq.octant_modifier = chi_mesh::Vector3(1.0, 1.0, 1.0);

      //==================================== Develop LDFE values
      DevelopSQLDFEValues(sq, legendre);

      initial_octant_SQs_.push_back(sq);
    } // for j
  }   // for i
}