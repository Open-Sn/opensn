#pragma once

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <filesystem>

/// Helper for building a MeshContinuum for an orthogonal mesh
/// given an array of nodes (one for each dimension)
std::shared_ptr<opensn::MeshContinuum>
BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets);

std::shared_ptr<opensn::MeshContinuum> BuildLineMesh(double length, unsigned int n, double xmin);

std::shared_ptr<opensn::MeshContinuum> BuildSquareMesh(double length, unsigned int n, double xmin);

std::shared_ptr<opensn::MeshContinuum> BuildBoxMesh(double length, unsigned int n, double xmin);

std::shared_ptr<opensn::MeshContinuum> BuildMeshFromFile(std::filesystem::path file_name);
