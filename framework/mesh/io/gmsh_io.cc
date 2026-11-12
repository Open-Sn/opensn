// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include <fstream>
#include <string>
#include <vector>

namespace opensn
{

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmsh(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromFile";

  // Open file
  std::ifstream file(options.file_name);
  if (not file.is_open())
    throw std::runtime_error(fname + ": Failed to open file " + options.file_name);

  // Check file format version
  std::string file_line;
  while (std::getline(file, file_line))
    if ("$MeshFormat" == file_line)
      break;

  std::getline(file, file_line);
  std::istringstream iss;
  iss = std::istringstream(file_line);

  file.close();

  double format = 0.0;
  if (not(iss >> format))
    throw std::logic_error(fname + ": Failed to read Gmsh file format.");
  else if (format == 2.2)
    return FromGmshV22(options);
  else if (format == 4.1)
    return FromGmshV41(options);
  else
    throw std::logic_error(fname + ": Only Gmsh formats 2.2 and 4.1 are supported.");

  return nullptr;
}

void
MeshIO::AddGmshPrismOrPyramidFaces(const int element_type, UnpartitionedMesh::LightWeightCell& cell)
{
  // Gmsh orders the nodes of both elements so that the first face (the triangle 0-1-2 of a prism,
  // the quadrilateral 0-1-2-3 of a pyramid) is counterclockwise seen from the opposite vertices,
  // as VTK does. The faces below are oriented outward for that ordering.
  const auto& v = cell.vertex_ids;
  std::vector<std::vector<uint64_t>> faces;
  if (element_type == 6)
  {
    OpenSnLogicalErrorIf(v.size() != 6, "MeshIO::FromGmsh: A prism needs 6 nodes.");
    faces = {{v[0], v[1], v[4], v[3]},
             {v[1], v[2], v[5], v[4]},
             {v[2], v[0], v[3], v[5]},
             {v[3], v[4], v[5]},
             {v[0], v[2], v[1]}};
  }
  else if (element_type == 7)
  {
    OpenSnLogicalErrorIf(v.size() != 5, "MeshIO::FromGmsh: A pyramid needs 5 nodes.");
    faces = {{v[0], v[3], v[2], v[1]},
             {v[0], v[1], v[4]},
             {v[1], v[2], v[4]},
             {v[2], v[3], v[4]},
             {v[3], v[0], v[4]}};
  }
  else
    OpenSnLogicalError("MeshIO::FromGmsh: Element type " + std::to_string(element_type) +
                       " is not a prism or pyramid.");
  for (auto& vertex_ids : faces)
  {
    UnpartitionedMesh::LightWeightFace face;
    face.vertex_ids = std::move(vertex_ids);
    cell.faces.push_back(std::move(face));
  }
}

} // namespace opensn
