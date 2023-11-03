#pragma once

#include <cstddef>
#include <vector>

namespace chi_mesh
{

//###################################################################
/**Utility class for handling face categorizations based on number
 * of vertices.*/
class GridFaceHistogram
{
private:
  std::vector<std::pair<size_t, size_t>> face_categories_;

public:
  explicit GridFaceHistogram(const std::vector<std::pair<size_t, size_t>>& face_categories)
    : face_categories_(face_categories)
  {
  }

  size_t NumberOfFaceHistogramBins() const;
  size_t MapFaceHistogramBins(size_t num_face_verts) const;
  size_t GetFaceHistogramBinDOFSize(size_t bin_number) const;
};

} // namespace chi_mesh


