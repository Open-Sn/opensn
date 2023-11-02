#include "chi_grid_face_histogram.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_mesh
{

/**Returns the number of bins.*/
size_t
GridFaceHistogram::NumberOfFaceHistogramBins() const
{
  return face_categories_.size();
}

/**Finds which bin holds the given number of vertices*/
size_t
GridFaceHistogram::MapFaceHistogramBins(size_t num_face_verts) const
{
  size_t category_counter = -1;
  for (auto category : face_categories_)
  {
    category_counter++;
    if (num_face_verts <= category.first) return category_counter;
  }

  return 0;
}

/**Finds the amount of vertices per face for the given bin.*/
size_t
GridFaceHistogram::GetFaceHistogramBinDOFSize(size_t bin_number) const
{
  size_t face_dof_size;

  try
  {
    face_dof_size = face_categories_.at(bin_number).first;
  }
  catch (std::out_of_range& o)
  {
    Chi::log.LogAllWarning() << "Fault detected in chi_mesh::MeshContinuum::"
                             << "GetFaceHistogramBinDOFSize.";
    return 0;
  }

  return face_dof_size;
}

} // namespace chi_mesh