#include "chi_surfacemesh.h"

#include "graphs/chi_directed_graph.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>

//###################################################################
/**Checks for cyclic dependencies in this mesh.
 * Transport type sweeps have a step where the inter-cell dependency
 * is acyclically sorted. This step is repeated here.*/
void
chi_mesh::SurfaceMesh::CheckCyclicDependencies(int num_angles)
{
  double tolerance = 1.0e-8;
  double dvarphi = 2.0 * M_PI / num_angles;
  chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  //================================================== Loop over angles
  for (int a = 0; a < num_angles; a++)
  {
    double varphi = 0.5 * dvarphi + a * dvarphi;

    chi_mesh::Vector3 omega;
    omega.x = cos(varphi);
    omega.y = sin(varphi);
    omega.z = 0.0;

    //================================= Add all polyfaces to graph
    chi::DirectedGraph G;
    size_t num_loc_cells = poly_faces_.size();
    for (size_t c = 0; c < num_loc_cells; c++)
      G.AddVertex();

    //================================= Now construct dependencies
    for (size_t c = 0; c < num_loc_cells; c++)
    {
      chi_mesh::PolyFace* face = poly_faces_[c];

      size_t num_edges = face->edges.size();
      for (size_t e = 0; e < num_edges; e++)
      {
        int v0i = face->edges[e][0];
        int v1i = face->edges[e][1];

        chi_mesh::Vector3 v01 = vertices_[v1i] - vertices_[v0i];
        chi_mesh::Vector3 n = v01.Cross(khat);
        n = n / n.Norm();

        double mu = omega.Dot(n);
        int neighbor = face->edges[e][2];
        if ((mu > (0.0 + tolerance)) and (neighbor >= 0))
        {
          //          boost::add_edge(c,neighbor,G);
          G.AddEdge(c, neighbor);
        } // if outgoing
      }   // for edge
    }     // for cell

    //================================================== Generic topological
    //                                                   sorting
    //    typedef boost::graph_traits<CHI_D_GRAPH>::vertex_descriptor gVertex;
    //
    //    boost::property_map<CHI_D_GRAPH, boost::vertex_index_t>::type
    //      index_map = get(boost::vertex_index, G);
    //
    //    std::vector<gVertex> sorted_list;
    //    try{
    //      boost::topological_sort(G,std::back_inserter(sorted_list));
    //    }
    //    catch (const boost::bad_graph& exc)
    //    {
    //      chi::log.LogAllError()
    //        << "Function CheckCyclicDependencies. Detected cyclic depency.";
    //     chi::Exit(EXIT_FAILURE);
    //    }

    auto topological_order = G.GenerateTopologicalSort();
    if (topological_order.empty())
    {
      Chi::log.LogAllError() << "Function CheckCyclicDependencies. Detected cyclic depency.";
      Chi::Exit(EXIT_FAILURE);
    }

    //================================= Cleanup
    G.Clear();
  } // for angles

  GetMeshStats();

  Chi::log.Log() << "Cyclic dependency check complete. No cycles or "
                 << "bad mesh elements were detected";
}

//###################################################################
/**Gets simple mesh statistics.
 */
void
chi_mesh::SurfaceMesh::GetMeshStats()
{
  std::vector<double> areas;
  std::vector<double> histo_bins;
  std::vector<int> histo;

  int num_negative_areas = 0;

  //============================================= Compute areas for
  //                                              each polyface
  size_t num_loc_cells = poly_faces_.size();
  areas.resize(num_loc_cells);
  double max_area = 0.0;
  for (size_t c = 0; c < num_loc_cells; c++)
  {
    chi_mesh::PolyFace* face = poly_faces_[c];

    size_t num_edges = face->edges.size();
    double area = 0.0;
    for (size_t e = 0; e < num_edges; e++)
    {
      int v0i = face->edges[e][0];
      int v1i = face->edges[e][1];

      chi_mesh::Vector3 v01 = vertices_[v1i] - vertices_[v0i];
      chi_mesh::Vector3 v02 = face->face_centroid - vertices_[v0i];

      // This is essentially the combine of the triangle for each side

      area += 0.5 * (v01.x * v02.y - v01.y * v02.x);
    } // for edge

    areas[c] = area;
    if (area > max_area) max_area = area;

    if (area <= 0.0) num_negative_areas += 1;
  } // for cell

  //============================================= Sort the areas
  std::sort(areas.begin(), areas.end(), std::greater<double>());

  //============================================= Compute histogram bins
  histo_bins.resize(10);
  histo.resize(10, 0);
  histo_bins[0] = max_area * 1.05;
  for (int i = 1; i < 10; i++)
    histo_bins[i] = histo_bins[i - 1] / 2.0;

  //============================================= Polulate histogram
  for (auto area : areas)
  {
    int home_bin = 9;
    for (int i = 0; i < 10; i++)
    {

      if (area <= histo_bins[i]) home_bin = i;

    } // check bins

    histo[home_bin] += 1;
  } // for areas

  std::stringstream output;
  for (int i = 0; i < 10; i++)
  {
    char buff[100];
    snprintf(buff, 100, "%11.3e", histo_bins[i]);

    output << "Areas < " << buff << " = " << histo[i] << "\n";
  }
  output << "Number of negative or zero faces = " << num_negative_areas;

  Chi::log.Log() << output.str();
}

//###################################################################
/**Computes load balancing parameters from a set of predictive cuts.
 * Does not actually perform these cuts.*/
void
chi_mesh::SurfaceMesh::ComputeLoadBalancing(std::vector<double>& x_cuts,
                                            std::vector<double>& y_cuts)
{
  Chi::log.Log() << "X-cuts to be logged: " << x_cuts.size();
  //  for (auto& val : x_cuts)
  //    chi::log.Log() << val;
  //
  Chi::log.Log() << "Y-cuts to be logged: " << y_cuts.size();
  //  for (auto& val : y_cuts)
  //    chi::log.Log() << val;

  //======================================== Sort faces into bins
  size_t I = x_cuts.size();
  size_t J = y_cuts.size();

  std::vector<std::vector<int>> IJ_bins(I + 1, std::vector<int>(J + 1, 0));

  for (auto& poly_face : poly_faces_)
  {
    int ref_i = 0;
    int ref_j = 0;
    for (size_t i = 0; i < I; ++i)
    {
      if (poly_face->face_centroid.x >= x_cuts[i]) ref_i = i + 1;
    } // for i
    for (size_t j = 0; j < J; ++j)
    {
      if (poly_face->face_centroid.y >= y_cuts[j]) ref_j = j + 1;
    } // for j

    IJ_bins[ref_i][ref_j] += 1;
  } // for face

  //======================================== Determine average and max
  int max_bin_size = 0;
  int tot_bin_size = 0;
  int i_max = 0, j_max = 0;

  for (int i = 0; i < (I + 1); ++i)
  {
    for (int j = 0; j < (J + 1); ++j)
    {
      if (IJ_bins[i][j] > max_bin_size)
      {
        max_bin_size = IJ_bins[i][j];
        i_max = i;
        j_max = j;
      }
      tot_bin_size += IJ_bins[i][j];
    }
  }

  double average = tot_bin_size / ((double)(I + 1) * (J + 1));

  Chi::log.Log() << "Average faces per set: " << average;
  Chi::log.Log() << "Maximum faces per set: " << max_bin_size << " at (i,j)= ( " << i_max << " , "
                 << j_max << " )";

  if (i_max == I) Chi::log.Log() << "X greater than " << x_cuts[i_max - 1];
  else if (i_max == 0)
    Chi::log.Log() << "X less than " << x_cuts[0];
  else
    Chi::log.Log() << "X greater than " << x_cuts[i_max - 1] << " and less than " << x_cuts[i_max];

  if (j_max == J) Chi::log.Log() << "Y greater than " << y_cuts[j_max - 1];
  else if (j_max == 0)
    Chi::log.Log() << "Y less than " << y_cuts[0];
  else
    Chi::log.Log() << "Y greater than " << y_cuts[j_max - 1] << " and less than " << y_cuts[j_max];

  Chi::log.Log() << "Max-to-average ratio: " << max_bin_size / average;
}
