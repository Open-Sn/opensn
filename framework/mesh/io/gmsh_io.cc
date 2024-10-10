// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <fstream>

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

  double format;
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

} // namespace opensn
