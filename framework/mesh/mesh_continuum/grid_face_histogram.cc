// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/grid_face_histogram.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

size_t
GridFaceHistogram::GetNumberOfFaceHistogramBins() const
{
  return face_categories_.size();
}

size_t
GridFaceHistogram::MapFaceHistogramBins(size_t num_face_verts) const
{
  int64_t category_counter = -1;
  for (auto category : face_categories_)
  {
    category_counter++;
    if (num_face_verts <= category.first)
      return category_counter;
  }

  return 0;
}

size_t
GridFaceHistogram::GetFaceHistogramBinDOFSize(size_t bin_number) const
{
  size_t face_dof_size = 0;

  try
  {
    face_dof_size = face_categories_.at(bin_number).first;
  }
  catch (std::out_of_range& o)
  {
    log.LogAllWarning() << "Fault detected in MeshContinuum::"
                        << "GetFaceHistogramBinDOFSize.";
    return 0;
  }

  return face_dof_size;
}

} // namespace opensn
