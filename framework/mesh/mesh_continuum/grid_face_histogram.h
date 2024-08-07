// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <vector>

namespace opensn
{

/// Utility class for handling face categorizations based on number of vertices.
class GridFaceHistogram
{
private:
  std::vector<std::pair<size_t, size_t>> face_categories_;

public:
  explicit GridFaceHistogram(const std::vector<std::pair<size_t, size_t>>& face_categories)
    : face_categories_(face_categories)
  {
  }

  /// Returns the number of bins.
  size_t NumberOfFaceHistogramBins() const;

  /// Finds which bin holds the given number of vertices
  size_t MapFaceHistogramBins(size_t num_face_verts) const;

  /// Finds the amount of vertices per face for the given bin.
  size_t GetFaceHistogramBinDOFSize(size_t bin_number) const;
};

} // namespace opensn
