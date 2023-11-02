#include "sldfe_sq.h"

//###################################################################
/**Computes the area of a cell. This routine uses
 * Girard's theorem to get the area of a spherical
 * triangle using the spherical excess.*/
double chi_math::SimplifiedLDFESQ::Quadrature::
  ComputeSphericalQuadrilateralArea(
    std::array<chi_mesh::Vertex,4>& vertices_xyz)
{
  const auto num_verts = 4;

  //======================================== Compute centroid
  chi_mesh::Vertex centroid_xyz;
  for (auto& v : vertices_xyz)
    centroid_xyz += v;
  centroid_xyz /= num_verts;

  //======================================== Compute area via
  //                                         spherical excess
  double area = 0.0;
  auto v0 = centroid_xyz.Normalized();
  for (int v=0; v<num_verts; ++v)
  {
    auto v1 = vertices_xyz[v];
    auto v2 = vertices_xyz[(v < (num_verts-1))? v+1 : 0];

    if ((v1-v0).Cross(v2-v0).Dot(v0) < 0.0)
      std::swap(v1,v2);

    //====================================== Lambda for spherical excess
    auto GetSphericalExcess = [](const chi_mesh::Vector3& vA,
                                 const chi_mesh::Vector3& vB,
                                 const chi_mesh::Vector3& vC)
    {
      const auto& n = vA;

      auto vAB = vB - vA;
      auto vAC = vC - vA;

      auto tAB = vAB.Cross(n).Normalized();
      auto tAC = vAC.Cross(n).Normalized();

      auto bAB = n.Cross(tAB).Normalized();
      auto bAC = n.Cross(tAC).Normalized();

      double mu = std::max(-1.0,std::fmin(1.0,bAB.Dot(bAC)));

      return std::fabs(acos(mu));
    };

    double excess = GetSphericalExcess(v0,v1,v2) +
                    GetSphericalExcess(v1,v2,v0) +
                    GetSphericalExcess(v2,v0,v1);

    area += excess - M_PI;
  }

  return area;
}


//###################################################################
/**Integrates shape functions to produce weights.*/
std::array<double,4> chi_math::SimplifiedLDFESQ::Quadrature::
IntegrateLDFEShapeFunctions(
  const SphericalQuadrilateral &sq,
  std::array<chi_math::DynamicVector<double>,4>& shape_coeffs,
  const std::vector<chi_math::QuadraturePointXYZ>& legendre_qpoints,
  const std::vector<double>& legendre_qweights)
{
  //=================================== Lambda to evaluate LDFE shape func
  auto EvaluateShapeFunction = [](chi_math::DynamicVector<double>& shape_coeffs,
                                  chi_mesh::Vector3& mu_eta_xi)
  {
    return shape_coeffs[0] + shape_coeffs[1]*mu_eta_xi[0] +
                             shape_coeffs[2]*mu_eta_xi[1] +
                             shape_coeffs[3]*mu_eta_xi[2];
  };

  //=================================== Determine integration bounds
  double x_tilde_max = 0.0;
  double x_tilde_min = 1.0;
  double y_tilde_max = 0.0;
  double y_tilde_min = 1.0;

  for (auto& v : sq.vertices_xy_tilde)
  {
    x_tilde_max = std::fmax(x_tilde_max,v.x);
    x_tilde_min = std::fmin(x_tilde_min,v.x);
    y_tilde_max = std::fmax(y_tilde_max,v.y);
    y_tilde_min = std::fmin(y_tilde_min,v.y);
  }

  //=================================== Integrate Legendre Quadrature
  std::array<double,4> integral = {0.0,0.0,0.0,0.0};
  int Nq = legendre_qpoints.size();
  double dx_tilde = (x_tilde_max - x_tilde_min);
  double dy_tilde = (y_tilde_max - y_tilde_min);

  for (int i=0; i<Nq; ++i)
  {
    for (int j=0; j<Nq; ++j)
    {
      //========================== Determine xy_tilde
      double x_tilde = x_tilde_min + (1.0 + legendre_qpoints[j][0])*dx_tilde/2.0;
      double y_tilde = y_tilde_min + (1.0 + legendre_qpoints[i][0])*dy_tilde/2.0;
      chi_mesh::Vector3 xy_tilde(x_tilde,y_tilde,0.0);

      //========================== Map to xyz
      auto xyz = (sq.rotation_matrix*xy_tilde + sq.translation_vector).Normalized();

      //========================== Determine Jacobian
      double r = sqrt(x_tilde*x_tilde + y_tilde*y_tilde + a*a);
      double detJ = (a/(r*r*r))*dx_tilde*dy_tilde/4.0;

      //========================== Evaluate shape funcs and add to integral
      for (int k=0; k<4; ++k)
        integral[k] += EvaluateShapeFunction(shape_coeffs[k], xyz)*detJ*
                       legendre_qweights[i]*legendre_qweights[j];
    }//for j
  }//for i

  return integral;
}

//###################################################################
/**Deploys the current set of SQs to all octants.*/
void chi_math::SimplifiedLDFESQ::Quadrature::CopyToAllOctants()
{
  deployed_SQs_.clear(); //just to be sure
  deployed_SQs_.reserve(initial_octant_SQs_.size() * 8);

  //======================================== Define modifying variables
  chi_mesh::Vector3 octant_mod(1.0,1.0,1.0);

  //======================================== Top NE octant, no change
  for (auto& sq : initial_octant_SQs_)
    deployed_SQs_.push_back(sq);

  //======================================== Top NW octant
  octant_mod = chi_mesh::Vector3(-1.0,1.0,1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Top SW octant
  octant_mod = chi_mesh::Vector3(-1.0,-1.0,1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Top SE octant
  octant_mod = chi_mesh::Vector3(1.0,-1.0,1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Bot NE octant
  octant_mod = chi_mesh::Vector3(1.0,1.0,-1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Bot NW octant
  octant_mod = chi_mesh::Vector3(-1.0,1.0,-1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Bot SW octant
  octant_mod = chi_mesh::Vector3(-1.0,-1.0,-1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Bot SE octant
  octant_mod = chi_mesh::Vector3(1.0,-1.0,-1.0);
  for (auto& sq : initial_octant_SQs_)
  {
    SphericalQuadrilateral new_sq = sq;

    for (auto& xyz : new_sq.vertices_xyz)   xyz = xyz*octant_mod;
    auto& vcc = new_sq.centroid_xyz;        vcc = vcc*octant_mod;
    for (auto& xyz : new_sq.sub_sqr_points) xyz = xyz*octant_mod;
    new_sq.octant_modifier                      = octant_mod;

    deployed_SQs_.push_back(new_sq);
  }

  //======================================== Make history entry
  deployed_SQs_history_.push_back(deployed_SQs_);
}

//###################################################################
/**Populates the quadrature abscissaes, weights and direction vectors.*/
void chi_math::SimplifiedLDFESQ::Quadrature::PopulateQuadratureAbscissae()
{
  abscissae_.clear();
  weights_.clear();
  omegas_.clear();

  for (const auto& sq : deployed_SQs_)
  {
    for (int i=0;i<4;++i)
    {
      const auto& omega = sq.sub_sqr_points[i];
      const double weight = sq.sub_sqr_weights[i];

      double theta = acos(omega.z);
      double phi = acos(omega.x/sin(theta));

      if (omega.y/sin(theta)<0.0)
        phi = 2.0*M_PI - phi;

      const auto abscissa = chi_math::QuadraturePointPhiTheta(phi, theta);

      abscissae_.push_back(abscissa);
      weights_.push_back(weight);
      omegas_.push_back(omega);
    }
  }
}

