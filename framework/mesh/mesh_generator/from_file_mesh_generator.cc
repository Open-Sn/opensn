// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/from_file_mesh_generator.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mesh/io/mesh_io.h"
#include "framework/utils/utils.h"
#include <filesystem>

namespace opensn
{

OpenSnRegisterObjectInNamespace(mesh, FromFileMeshGenerator);

InputParameters
FromFileMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription("Generator for loading an unpartitioned mesh from a file.");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddRequiredParameter<std::string>("filename", "Path to the file.");
  params.AddOptionalParameter("material_id_fieldname",
                              "BlockID",
                              "The name of the field storing cell block/material ids. Only really "
                              "used for .vtu, .pvtu and .e files.");
  params.AddOptionalParameter(
    "boundary_id_fieldname", "", "The name of the field storing boundary-ids");

  return params;
}

FromFileMeshGenerator::FromFileMeshGenerator(const InputParameters& params)
  : MeshGenerator(params),
    filename_(params.GetParamValue<std::string>("filename")),
    material_id_fieldname_(params.GetParamValue<std::string>("material_id_fieldname")),
    boundary_id_fieldname_(params.GetParamValue<std::string>("boundary_id_fieldname"))
{
}

std::shared_ptr<UnpartitionedMesh>
FromFileMeshGenerator::GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh)
{
  OpenSnInvalidArgumentIf(input_umesh != nullptr,
                          "FromFileMeshGenerator can not be preceded by another"
                          " mesh generator because it cannot process an input mesh");

  UnpartitionedMesh::Options options;
  options.file_name = filename_;
  options.scale = scale_;
  options.material_id_fieldname = material_id_fieldname_;
  options.boundary_id_fieldname = boundary_id_fieldname_;

  const std::filesystem::path filepath(filename_);
  AssertReadableFile(filename_);

  log.Log() << "FromFileMeshGenerator: Generating UnpartitionedMesh";
  const std::string extension = filepath.extension();
  if (extension == ".obj")
    return MeshIO::FromOBJ(options);
  else if (extension == ".msh")
    return MeshIO::FromGmsh(options);
  else if (extension == ".e")
    return MeshIO::FromExodus(options);
  else if (extension == ".vtu")
    return MeshIO::FromVTU(options);
  else if (extension == ".pvtu")
    return MeshIO::FromPVTU(options);
  else if (extension == ".case")
    return MeshIO::FromEnsightGold(options);
  else
    OpenSnInvalidArgument("Unsupported file type \"" + extension +
                          "\". Supported types limited to"
                          ".obj, .msh, .e, .vtu, .pvtu, .case.");
}

} // namespace opensn
