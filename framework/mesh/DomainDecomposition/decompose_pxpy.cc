#include "opensn/framework/mesh/SurfaceMesh/chi_surfacemesh.h"

#include <algorithm>

typedef std::pair<int, int> IntPair;
typedef std::vector<double> DblVec;
typedef std::vector<int> IntVec;

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

//================================================== Define LBF-Calc funtion
/**Makes a centroid based load balance factor calculation.
 *
 * \author Jan*/
double
chi_mesh::ComputeLBF(std::vector<Vector3>& points,
                     std::vector<double>& x_cuts,
                     std::vector<double>& y_cuts)
{
  std::vector<IntVec> balance_matrix(x_cuts.size() + 1, IntVec(y_cuts.size() + 1, 0));

  size_t num_points = points.size();
  size_t num_x_cuts = x_cuts.size();
  size_t num_y_cuts = y_cuts.size();
  for (size_t p = 0; p < num_points; p++)
  {
    double xr = points[p].x;
    double yr = points[p].y;

    int ir = -1;
    for (int i = 0; i < num_x_cuts; i++)
    {
      if (xr <= x_cuts[i])
      {
        ir = i;
        break;
      }
    }
    if (ir < 0) ir = num_x_cuts;

    int jr = -1;
    for (int j = 0; j < num_y_cuts; j++)
    {
      if (yr <= y_cuts[j])
      {
        jr = j;
        break;
      }
    }
    if (jr < 0) jr = num_y_cuts;

    balance_matrix[ir][jr] += 1;
  }

  // Now find max and average
  int max_points_per_bin = 0;
  int bin_sum = 0;
  for (int i = 0; i <= num_x_cuts; i++)
  {
    for (int j = 0; j <= num_y_cuts; j++)
    {
      bin_sum += balance_matrix[i][j];

      if (balance_matrix[i][j] > max_points_per_bin) max_points_per_bin = balance_matrix[i][j];
    } // for j
  }   // for i
  double average = bin_sum / (num_x_cuts + 1) / (num_y_cuts + 1);

  return max_points_per_bin / average;
}

//###################################################################
/** Decomposes a 2D surface mesh using the centroids in a Px-Py fashion.*/
void
chi_mesh::DecomposeSurfaceMeshPxPy(const chi_mesh::SurfaceMesh& smesh, int px, int py)
{
  //================================================== Collect centroids
  size_t num_pfaces = smesh.GetPolygons().size();
  std::vector<chi_mesh::Vector3> centroids(num_pfaces);
  for (int pf = 0; pf < num_pfaces; pf++)
    centroids[pf] = smesh.GetPolygons()[pf]->face_centroid;

  //================================================== Define sort operators
  struct
  {
    bool operator()(chi_mesh::Vector3 a, chi_mesh::Vector3 b)
    {
      if (a.x < b.x) { return true; }
      else { return false; }
    }
  } compare_x;

  struct
  {
    bool operator()(chi_mesh::Vector3 a, chi_mesh::Vector3 b)
    {
      if (a.y < b.y) { return true; }
      else { return false; }
    }
  } compare_y;

  //================================================== Create sorts
  std::vector<chi_mesh::Vector3> centroids_sortedx;
  std::vector<chi_mesh::Vector3> centroids_sortedy;
  std::copy(centroids.begin(), centroids.end(), std::back_inserter(centroids_sortedx));

  std::stable_sort(centroids_sortedx.begin(), centroids_sortedx.end(), compare_x);

  std::vector<IntPair> x_bins(px);
  std::vector<double> x_cuts(px - 1);
  std::vector<DblVec> y_cuts_per_x_bin(px, DblVec(py - 1));

  //================================================== Populate xbins
  int dx = std::ceil(centroids.size() / (double)px);
  for (int x = 0; x < px; x++)
  {
    x_bins[x].first = x * dx;
    x_bins[x].second = x * dx + dx - 1;

    if (x == (px - 1)) x_bins[x].second = centroids_sortedx.size() - 1;
  }

  //================================================== Populate x-cuts
  for (int x = 0; x < (px - 1); x++)
  {
    int up_bounds_x = x * dx + dx - 1;
    double a = centroids_sortedx[up_bounds_x].x;
    double b = centroids_sortedx[up_bounds_x + 1].x;

    x_cuts[x] = (b + a) / 2.0;

    Chi::log.Log() << "X-cut" << x << " " << x_cuts[x];
  }

  //================================================== Balance according to each
  //                                                   x-bin
  double min_lbf = 1000.0; // Minimum load balance factor
  int min_bin = 0;
  for (int x = 0; x < px; x++)
  {
    int bin_i = x_bins[x].first;
    int bin_f = x_bins[x].second;

    centroids_sortedy.clear();
    for (int i = bin_i; i <= bin_f; i++)
      centroids_sortedy.push_back(centroids_sortedx[i]);

    std::stable_sort(centroids_sortedy.begin(), centroids_sortedy.end(), compare_y);

    int dy = std::ceil(centroids_sortedy.size() / (double)py);
    for (int y = 0; y < (py - 1); y++)
    {
      int up_bounds_y = y * dy + dy - 1;

      double a = centroids_sortedy[up_bounds_y].y;
      double b = centroids_sortedy[up_bounds_y + 1].y;

      y_cuts_per_x_bin[x][y] = (b + a) / 2.0;
    }
    double lbf = ComputeLBF(centroids, x_cuts, y_cuts_per_x_bin[x]);

    if (lbf < min_lbf)
    {
      min_lbf = lbf;
      min_bin = x;
    }
    Chi::log.Log() << "Load balance factor: " << lbf;
  } // for x

  //================================================== Write y-cuts
  for (int y = 0; y < (py - 1); y++)
  {
    Chi::log.Log() << "Y-cut" << y << " " << y_cuts_per_x_bin[min_bin][y];
  }
}
