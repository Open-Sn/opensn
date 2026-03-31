#pragma once

#include "framework/mesh/mesh_continuum/mesh_continuum.h"

/// Helper for building a MeshContinuum for an orthogonal mesh
/// given an array of nodes (one for each dimension)
std::shared_ptr<opensn::MeshContinuum>
BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets);

std::shared_ptr<opensn::MeshContinuum> BuildSquareMesh(double length, unsigned int n, double xmin);
