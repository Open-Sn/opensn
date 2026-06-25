#pragma once

#include "framework/mesh/mesh/mesh.h"
#include <filesystem>

/// Helper for building a Mesh for an orthogonal mesh
/// given an array of nodes (one for each dimension)
std::shared_ptr<opensn::Mesh>
BuildOrthogonalMesh(const std::vector<std::vector<double>>& node_sets);

std::shared_ptr<opensn::Mesh> BuildLineMesh(double length, unsigned int n, double xmin);

std::shared_ptr<opensn::Mesh> BuildSquareMesh(double length, unsigned int n, double xmin);

std::shared_ptr<opensn::Mesh> BuildBoxMesh(double length, unsigned int n, double xmin);

std::shared_ptr<opensn::Mesh> BuildMeshFromFile(std::filesystem::path file_name);
