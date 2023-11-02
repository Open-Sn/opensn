#include "sldfe_sq.h"

#include <map>

//###################################################################
/**Develops LDFE quantities.*/
void chi_math::SimplifiedLDFESQ::Quadrature::
  DevelopSQLDFEValues(SphericalQuadrilateral& sq,
                      chi_math::QuadratureGaussLegendre& legendre)
{
  //============================================= Determine sq tilde center
  chi_mesh::Vertex sq_tilde_center;
  for (const auto& v : sq.vertices_xy_tilde)
    sq_tilde_center += v;
  sq_tilde_center/=4;

  //============================================= Determine off-set vectors
  auto& vc = sq_tilde_center;
  std::array<chi_mesh::Vector3,4> vctoi;
  for (int v=0;v<4;++v)
    vctoi[v] = sq.vertices_xy_tilde[v] - vc;

  //============================================= Determine sub-sub-squares
  std::array<std::array<chi_mesh::Vertex,4>,4> sub_sub_square_xy_tilde;
  std::map<std::string,chi_mesh::Vertex> vm;

  for (int v=0;v<4;++v) vm[std::to_string(v)] = sq.vertices_xy_tilde[v];

  vm["01"] = 0.5*(sq.vertices_xy_tilde[0]+sq.vertices_xy_tilde[1]);
  vm["12"] = 0.5*(sq.vertices_xy_tilde[1]+sq.vertices_xy_tilde[2]);
  vm["23"] = 0.5*(sq.vertices_xy_tilde[2]+sq.vertices_xy_tilde[3]);
  vm["03"] = 0.5*(sq.vertices_xy_tilde[0]+sq.vertices_xy_tilde[3]);
  vm["c"]  = sq_tilde_center;

  auto& sst = sub_sub_square_xy_tilde;
  sst[0] = {vm["0"],vm["01"],vm["c"],vm["03"]};
  sst[1] = {vm["01"],vm["1"],vm["12"],vm["c"]};
  sst[2] = {vm["c"],vm["12"],vm["2"],vm["23"]};
  sst[3] = {vm["03"],vm["c"],vm["23"],vm["3"]};

  //============================================= Determine sub-sub-square
  //                                              xyz
  std::array<std::array<chi_mesh::Vertex,4>,4> sub_sub_square_xyz;
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      sub_sub_square_xyz[i][j] =
        (sq.rotation_matrix*sub_sub_square_xy_tilde[i][j] +
         sq.translation_vector).Normalized();

  //============================================= Compute sub-sub-square area
  std::array<double,4> SA_i = {0.0,0.0,0.0,0.0};
  for (int i=0; i<4; ++i)
    SA_i[i] = ComputeSphericalQuadrilateralArea(sub_sub_square_xyz[i]);

  //============================================= Apply optimization
  if (qp_optimization_type_ == QuadraturePointOptimization::CENTROID)
    for (int i=0; i<4; ++i)
    {
      for (int j=0; j<4; ++j)
        sq.sub_sqr_points[i] += sub_sub_square_xyz[i][j];
      sq.sub_sqr_points[i] /= 4.0;
      sq.sub_sqr_points[i].Normalize();

      sq.sub_sqr_weights[i] = SA_i[i];
    }//for i
  else if (qp_optimization_type_ == QuadraturePointOptimization::EMPIRICAL)
    EmpiricalQPOptimization(sq,legendre,vc,vctoi,SA_i);
  else if (qp_optimization_type_ == QuadraturePointOptimization::ISOLATED)
    IsolatedQPOptimization(sq,legendre,vc,vctoi,SA_i);

}

