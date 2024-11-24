// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

/// Structure containing edge properties
struct Edge
{
  /// Indices of the vertices
  std::array<int, 2> v_index{-1, -1};
  /// Indices of faces adjoining it
  std::array<int, 4> f_index{-1, -1, -1, -1};
  /// Vector vertices
  std::array<Vector3, 2> vertices{};
};

} // namespace opensn
