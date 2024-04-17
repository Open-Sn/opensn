// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"

namespace opensn
{

class GmshIO
{
public:
  static std::shared_ptr<UnpartitionedMesh> FromFile(const UnpartitionedMesh::Options& options);
};

} // namespace opensn
