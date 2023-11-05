#include "opensn/framework/mesh/MeshGenerator/FromFileMeshGenerator.h"

#include "opensn/framework/mesh/UnpartitionedMesh/unpartitioned_mesh.h"

#include "opensn/framework/ChiObjectFactory.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

#include <filesystem>

namespace chi_mesh
{

RegisterChiObject(chi_mesh, FromFileMeshGenerator);

chi::InputParameters
FromFileMeshGenerator::GetInputParameters()
{
  chi::InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription("Generator for loading an unpartitioned mesh"
                               " from a file.");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddRequiredParameter<std::string>("filename", "Path to the file.");
  params.AddOptionalParameter(
    "material_id_fieldname",
    "BlockID",
    "The name of the field storing cell block/material ids. Only really used "
    "for .vtu, .pvtu and .e files.");
  params.AddOptionalParameter(
    "boundary_id_fieldname", "", "The name of the field storing boundary-ids");

  return params;
}

FromFileMeshGenerator::FromFileMeshGenerator(const chi::InputParameters& params)
  : MeshGenerator(params),
    filename_(params.GetParamValue<std::string>("filename")),
    material_id_fieldname_(params.GetParamValue<std::string>("material_id_fieldname")),
    boundary_id_fieldname_(params.GetParamValue<std::string>("boundary_id_fieldname"))
{
  const std::filesystem::path filepath(filename_);
  const std::string extension = filepath.extension();
}

std::unique_ptr<UnpartitionedMesh>
FromFileMeshGenerator::GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh)
{
  ChiInvalidArgumentIf(input_umesh != nullptr,
                       "FromFileMeshGenerator can not be preceded by another"
                       " mesh generator because it cannot process an input mesh");

  UnpartitionedMesh::Options options;
  options.file_name = filename_;
  options.scale = scale_;
  options.material_id_fieldname = material_id_fieldname_;
  options.boundary_id_fieldname = boundary_id_fieldname_;

  const std::filesystem::path filepath(filename_);
  const std::string extension = filepath.extension();

  auto umesh = std::make_unique<UnpartitionedMesh>();

  Chi::log.Log() << "FromFileMeshGenerator: Generating UnpartitionedMesh";

  if (extension == ".obj") umesh->ReadFromWavefrontOBJ(options);
  else if (extension == ".msh")
    umesh->ReadFromMsh(options);
  else if (extension == ".e")
    umesh->ReadFromExodus(options);
  else if (extension == ".vtu")
    umesh->ReadFromVTU(options);
  else if (extension == ".pvtu")
    umesh->ReadFromPVTU(options);
  else if (extension == ".case")
    umesh->ReadFromEnsightGold(options);
  else
    ChiInvalidArgument("Unsupported file type \"" + extension +
                       "\". Supported types limited to"
                       ".obj, .msh, .e, .vtu, .pvtu, .case.");

  Chi::log.Log() << "FromFileMeshGenerator: Done generating UnpartitionedMesh";
  return umesh;
}

} // namespace chi_mesh
