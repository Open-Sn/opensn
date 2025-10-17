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

FromFileMeshGenerator::FromFileMeshGenerator(const InputParameters& params)
  : MeshGenerator(params),
    filename_(params.GetParamValue<std::string>("filename")),
    block_id_fieldname_(params.GetParamValue<std::string>("block_id_fieldname")),
    boundary_id_fieldname_(params.GetParamValue<std::string>("boundary_id_fieldname")),
    coord_sys_(params.GetParamValue<std::string>("coord_sys") == "cartesian"     ? CARTESIAN
               : params.GetParamValue<std::string>("coord_sys") == "cylindrical" ? CYLINDRICAL
                                                                                 : SPHERICAL)
{
}

std::shared_ptr<UnpartitionedMesh>
FromFileMeshGenerator::GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh)
{
  if (input_umesh != nullptr)
    throw std::invalid_argument("FromFileMeshGenerator can not be preceded by another "
                                "mesh generator because it cannot process an input mesh");

  UnpartitionedMesh::Options options;
  options.file_name = filename_;
  options.scale = scale_;
  options.block_id_fieldname = block_id_fieldname_;
  options.boundary_id_fieldname = boundary_id_fieldname_;

  const std::filesystem::path filepath(filename_);
  AssertReadableFile(filename_);

  log.Log() << "FromFileMeshGenerator: Generating UnpartitionedMesh";
  const std::string extension = filepath.extension();

  std::shared_ptr<UnpartitionedMesh> umesh;
  if (extension == ".obj")
    umesh = MeshIO::FromOBJ(options);
  else if (extension == ".msh")
    umesh = MeshIO::FromGmsh(options);
  else if (extension == ".e")
    umesh = MeshIO::FromExodusII(options);
  else if (extension == ".vtu")
    umesh = MeshIO::FromVTU(options);
  else if (extension == ".pvtu")
    umesh = MeshIO::FromPVTU(options);
  else if (extension == ".case")
    umesh = MeshIO::FromEnsightGold(options);
  else if (std::filesystem::is_directory(filepath)
           // Contains boundary patch information.
           && std::filesystem::exists(filepath / "constant/polyMesh/boundary")
           // Contains information for the general mesh
           && std::filesystem::exists(filepath / "constant/polyMesh/faces") &&
           std::filesystem::exists(filepath / "constant/polyMesh/neighbour") &&
           std::filesystem::exists(filepath / "constant/polyMesh/owner") &&
           std::filesystem::exists(filepath / "constant/polyMesh/points"))
    umesh = MeshIO::FromOpenFOAM(options);
  else
    throw std::invalid_argument("Unsupported file type \"" + extension +
                                "\". Supported types limited to "
                                ".obj, .msh, .e, .vtu, .pvtu, .case.");

  umesh->SetCoordinateSystem(coord_sys_);
  return umesh;
}

OpenSnRegisterObjectInNamespace(mesh, FromFileMeshGenerator);

InputParameters
FromFileMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription("Generator for loading an unpartitioned mesh from a file.");

  params.AddRequiredParameter<std::string>("filename", "Path to the file.");
  params.AddOptionalParameter("block_id_fieldname",
                              "BlockID",
                              "The name of the field storing cell block/material ids. Only really "
                              "used for .vtu, .pvtu and .e files.");
  params.AddOptionalParameter(
    "boundary_id_fieldname", "", "The name of the field storing boundary-ids");

  params.AddOptionalParameter("coord_sys", "cartesian", "The coordinate system of the mesh.");
  params.ConstrainParameterRange(
    "coord_sys", AllowableRangeList::New({"cartesian", "cylindrical", "spherical"}));

  return params;
}

std::shared_ptr<FromFileMeshGenerator>
FromFileMeshGenerator::Create(const ParameterBlock& params)
{
  const auto& factory = ObjectFactory::GetInstance();
  auto ptr = factory.Create<FromFileMeshGenerator>("mesh::FromFileMeshGenerator", params);
  return ptr;
}

} // namespace opensn
