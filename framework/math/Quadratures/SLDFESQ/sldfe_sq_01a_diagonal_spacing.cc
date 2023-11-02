#include "sldfe_sq.h"

//###################################################################
/**Generates diagonal spacings.*/
void
chi_math::SimplifiedLDFESQ::Quadrature::GenerateDiagonalSpacings(int level)
{
  //======================================== Define constants
  const int Ns = (level + 1); // Number of subdivisions
  const int Np = Ns + 1;      // Number of diagonal points

  const auto ihat = chi_mesh::Vector3(1.0, 0.0, 0.0);

  //======================================== Define rotation matrix
  chi_mesh::Matrix3x3 Rihat;
  auto n = chi_mesh::Vector3(0.0, -1.0 / sqrt(2), 1.0 / sqrt(2));
  auto& t = ihat;
  auto b = n.Cross(t).Normalized();

  Rihat.SetColJVec(0, t);
  Rihat.SetColJVec(1, b);
  Rihat.SetColJVec(2, n);

  //======================================== Generate p-points
  std::vector<chi_mesh::Vector3> p_points(Np);
  double dphi = acos(a) / Ns;
  double alpha = 0.10005;
  double beta = 1.0185;

  for (int i = 0; i < Np; ++i)
  {
    double phi = i * dphi * (1.0 + alpha * (cos(beta * M_PI_2 * i / Ns) - cos(beta * M_PI_2)));

    p_points[i] = Rihat * chi_mesh::Vector3(cos(phi), sin(phi), 0.0);
  }

  //======================================== Compute tilde points
  std::vector<chi_mesh::Vector3> tilde_points(Np);
  for (int i = 0; i < Np; ++i)
  {
    double R = a / p_points[i].x;
    double x_tilde = p_points[i].y * R;
    double y_tilde = p_points[i].z * R;

    tilde_points[i] = chi_mesh::Vector3(x_tilde, y_tilde, 0.0);
  }

  diagonal_vertices_ = tilde_points;
}