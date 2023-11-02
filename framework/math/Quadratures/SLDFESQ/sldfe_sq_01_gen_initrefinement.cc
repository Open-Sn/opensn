#include "sldfe_sq.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

//###################################################################
/**Generates uniform spherical quadrilaterals from the
 * subdivision of an inscribed cube.*/
void chi_math::SimplifiedLDFESQ::Quadrature::GenerateInitialRefinement(int level)
{
  chi::Timer timer;
  timer.Reset();
  initial_level_ = level;

  //======================================== Define constants
  const chi_mesh::Vector3 ihat = chi_mesh::Vector3(1.0,0.0,0.0);
  const chi_mesh::Vector3 jhat = chi_mesh::Vector3(0.0,1.0,0.0);
  const chi_mesh::Vector3 khat = chi_mesh::Vector3(0.0,0.0,1.0);

  //======================================== Build rotation matrices
  //                                         for cube faces
  chi_mesh::Matrix3x3 Rxface;
  Rxface.SetColJVec(0,jhat);
  Rxface.SetColJVec(1,khat);
  Rxface.SetColJVec(2,ihat);

  chi_mesh::Matrix3x3 Ryface;
  Ryface.SetColJVec(0,ihat);
  Ryface.SetColJVec(1,khat);
  Ryface.SetColJVec(2,-1.0*jhat);

  chi_mesh::Matrix3x3 Rzface;
  Rzface.SetDiagonalVec(1.0,1.0,1.0);

  //======================================== Set translation vectors
  //                                         for cube faces
  auto txface = a*ihat;
  auto tyface = a*jhat;
  auto tzface = a*khat;

  //======================================== Generate general diagonal
  //                                         spacings in xy-tilde coordinates
  GenerateDiagonalSpacings(level);

  //======================================== Generate vertices for each face
  //                                         of inscribed cube
  GenerateReferenceFaceVertices(Rxface,txface,level);
  GenerateReferenceFaceVertices(Ryface,tyface,level);
  GenerateReferenceFaceVertices(Rzface,tzface,level);

  //======================================== Compute areas
  double total_area = 0.0;
  double area_max = -100.0;
  double area_min =  100.0;
  bool negative_weights_found = false;
  for (auto& sq : initial_octant_SQs_)
  {
    double area = 0.0;
    for (int i=0; i<4; ++i)
    {
      area += sq.sub_sqr_weights[i];
      if (area < 0.0) negative_weights_found = true;
    }
    total_area += area;
    area_max = std::fmax(area_max,area);
    area_min = std::fmin(area_min,area);
  }
  double area_avg = total_area / initial_octant_SQs_.size();

  if (negative_weights_found)
    Chi::log.Log0Warning()
      << "SLDFESQ Quadrature detected negative weights.";

  //======================================== Print Statistics
  double time = timer.GetTime()/1000.0;
  Chi::log.Log0Verbose1() << "Number of dirs/octant: " << initial_octant_SQs_.size();
  Chi::log.Log0Verbose1() << "Total weight         : " << total_area;
  Chi::log.Log0Verbose1() << "Total weight/(pi/2)  : " << total_area/M_PI_2;
  Chi::log.Log0Verbose1() << "Area Max/Min         : " << area_max/area_min;
  Chi::log.Log0Verbose1() << "Area Max/Avg         : " << area_max/area_avg;

  CopyToAllOctants();

  //======================================== Populate quadriture points
  PopulateQuadratureAbscissae();

  Chi::log.Log0Verbose1() << "Time taken           : " << time;
}