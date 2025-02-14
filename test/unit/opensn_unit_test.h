#pragma once

#include <gtest/gtest.h>
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

using namespace opensn;

class OpenSnUnitTest : public ::testing::Test
{
public:
  void SetUp() override;
  void TearDown() override;

protected:
  /// Helper for building a MeshContinuum for an orthogonal mesh
  /// given an array of nodes (one for each dimension)
  std::shared_ptr<MeshContinuum>
  BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets) const;
};
