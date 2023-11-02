#include "sldfe_sq.h"

#include <map>
#include "chi_runtime.h"
#include "chi_log.h"

#include "mesh/chi_meshvector.h"

//###################################################################
/**Split a SQ.*/
std::array<chi_math::SimplifiedLDFESQ::SphericalQuadrilateral,4>
  chi_math::SimplifiedLDFESQ::Quadrature::
    SplitSQ(SphericalQuadrilateral &sq,
            chi_math::QuadratureGaussLegendre& legendre)
{
  std::array<SphericalQuadrilateral,4> new_sqs;

  //============================================= Determine sq tilde center
  chi_mesh::Vertex sq_tilde_center;
  for (const auto& v : sq.vertices_xy_tilde)
    sq_tilde_center += v;
  sq_tilde_center/=4;

  //============================================= Determine sub-sub-squares
  //                                              Tilde coordinates
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

  for (int i=0; i<4; ++i)
    new_sqs[i].vertices_xy_tilde = sst[i];

  //============================================= Determine xyz-prime
  for (int i=0; i<4; ++i)
    for (int v=0; v<4; ++v)
      new_sqs[i].vertices_xyz_prime[v] = sq.rotation_matrix*sst[i][v] +
                                         sq.translation_vector;

  //============================================= Compute xyz
  for (int i=0; i<4; ++i)
    for (int v=0; v<4; ++v)
      new_sqs[i].vertices_xyz[v] = new_sqs[i].vertices_xyz_prime[v].Normalized();

  //============================================= Compute SQ xyz-centroid,
  //                                              R,T,area, ldfe
  for (int i=0; i<4; ++i)
  {
    for (int v=0; v<4; ++v)
      new_sqs[i].centroid_xyz += new_sqs[i].vertices_xyz[v];
    new_sqs[i].centroid_xyz /= 4;
    new_sqs[i].centroid_xyz = new_sqs[i].centroid_xyz.Normalized()*sq.octant_modifier;

    new_sqs[i].rotation_matrix = sq.rotation_matrix;
    new_sqs[i].translation_vector = sq.translation_vector;

    new_sqs[i].area = ComputeSphericalQuadrilateralArea(new_sqs[i].vertices_xyz);
    DevelopSQLDFEValues(new_sqs[i],legendre);
    new_sqs[i].octant_modifier = sq.octant_modifier;

    for (int v=0; v<4; ++v)
    {
      new_sqs[i].vertices_xyz[v] = new_sqs[i].vertices_xyz[v]*sq.octant_modifier;
      new_sqs[i].sub_sqr_points[v] = new_sqs[i].sub_sqr_points[v]*sq.octant_modifier;
    }
  }

  return new_sqs;
}

//###################################################################
/**Locally refines the cells.*/
void chi_math::SimplifiedLDFESQ::Quadrature::
  LocallyRefine(const chi_mesh::Vector3 &ref_dir,
                const double cone_size,
                const bool dir_as_plane_normal)
{
  auto ref_dir_n = ref_dir.Normalized();
  double mu_cone = cos(cone_size);
  std::vector<SphericalQuadrilateral> new_deployment;
  new_deployment.reserve(deployed_SQs_.size());

  chi_math::QuadratureGaussLegendre legendre(QuadratureOrder::THIRTYSECOND);

  int num_refined = 0;
  for (auto& sq : deployed_SQs_)
  {
    bool sq_to_be_split = false;

    if (not dir_as_plane_normal)
      sq_to_be_split = sq.centroid_xyz.Dot(ref_dir_n)>mu_cone;
    else
      sq_to_be_split = std::fabs(sq.centroid_xyz.Dot(ref_dir_n)) <
                       (sin(cone_size));

    if (not sq_to_be_split)
      new_deployment.push_back(sq);
    else
    {
      auto new_sqs = SplitSQ(sq,legendre);
      for (auto& nsq : new_sqs)
        new_deployment.push_back(nsq);
      ++num_refined;
    }
  }

  deployed_SQs_.clear();
  deployed_SQs_ = new_deployment;
  deployed_SQs_history_.push_back(new_deployment);

  PopulateQuadratureAbscissae();

  Chi::log.Log() << "SLDFESQ refined " << num_refined << " SQs.";
}