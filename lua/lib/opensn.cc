// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/logging/log.h"
#include "framework/materials/material_property.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/mesh/mesh_vector.h"
#include "framework/parameters/input_parameters.h"
#include "framework/parameters/parameter_block.h"
#include "framework/post_processors/post_processor.h"
#include "lua/lib/console.h"
#include "lua/lib/types.h"
#include "lua/lib/enum.h"
#include "lua/lib/functions.h"
#include "framework/runtime.h"
#include "framework/math/quadratures/gauss_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/gausschebyshev_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/functions/function.h"
#include "framework/math/functions/scalar_spatial_function.h"
#include "framework/materials/material.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/isotropic_multigroup_source.h"
#include "framework/mesh/logical_volume/boolean_logical_volume.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/logical_volume/rcc_logical_volume.h"
#include "framework/mesh/logical_volume/rpp_logical_volume.h"
#include "framework/mesh/logical_volume/sphere_logical_volume.h"
#include "framework/mesh/logical_volume/surface_mesh_logical_volume.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"
#include "framework/mesh/mesh_generator/from_file_mesh_generator.h"
#include "framework/mesh/mesh_generator/extruder_mesh_generator.h"
#include "framework/mesh/mesh_generator/split_file_mesh_generator.h"
#include "framework/mesh/mesh_generator/distributed_mesh_generator.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/graphs/linear_graph_partitioner.h"
#include "framework/graphs/petsc_graph_partitioner.h"
#include "framework/field_functions/field_function.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "framework/physics/solver.h"
#include "framework/post_processors/aggregate_nodal_value_post_processor.h"
#include "framework/post_processors/cell_volume_integral_post_processor.h"
#include "framework/post_processors/solver_info_post_processor.h"
#include "modules/diffusion/diffusion_solver.h"
#include "modules/diffusion/cfem_diffusion_solver.h"
#include "modules/diffusion/dfem_diffusion_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_solver/lbs_curvilinear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/lbs_mip_solver.h"
#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"
#include "modules/linear_boltzmann_solvers/executors/nl_keigen.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_scdsa.h"
#include "modules/linear_boltzmann_solvers/executors/pi_keigen_smm.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/response_evaluator.h"
#include "modules/point_reactor_kinetics/point_reactor_kinetics.h"
#include "lua/lib/aquad.h"
#include "lua/lib/mesh.h"
#include "lua/lib/fieldfunc.h"
#include "lua/lib/logger.h"
#include "lua/lib/post.h"
#include "lua/lib/solver.h"
#include "lua/lib/xs.h"
#include "lua/lib/mat.h"
#include "lua/lib/lbs.h"
#include "LuaBridge/LuaBridge.h"
#include <lua.h>
#include <memory>
#include <stdexcept>

using namespace opensn;

namespace opensnlua
{

void
MPIBarrier()
{
  opensn::mpi_comm.barrier();
}

static bool reg = opensnlua::Console::Bind(
  [](lua_State* L)
  {
    luabridge::getGlobalNamespace(L)
      .addFunction("MPIBarrier", &opensnlua::MPIBarrier)
      .addFunction("Exit", &Exit);

    // enum values
    luabridge::getGlobalNamespace(L)
      .addVariable("LOG_0", Logger::LOG_0)
      .addVariable("LOG_0WARNING", Logger::LOG_0WARNING)
      .addVariable("LOG_0ERROR", Logger::LOG_0ERROR)
      .addVariable("LOG_0VERBOSE_0", Logger::LOG_0VERBOSE_0)
      .addVariable("LOG_0VERBOSE_1", Logger::LOG_0VERBOSE_1)
      .addVariable("LOG_0VERBOSE_2", Logger::LOG_0VERBOSE_2)
      .addVariable("LOG_ALL", Logger::LOG_ALL)
      .addVariable("LOG_ALLWARNING", Logger::LOG_ALLWARNING)
      .addVariable("LOG_ALLERROR", Logger::LOG_ALLERROR)
      .addVariable("LOG_ALLVERBOSE_0", Logger::LOG_ALLVERBOSE_0)
      .addVariable("LOG_ALLVERBOSE_1", Logger::LOG_ALLVERBOSE_1)
      .addVariable("LOG_ALLVERBOSE_2", Logger::LOG_ALLVERBOSE_2);

    luabridge::getGlobalNamespace(L)
      .addVariable("OP_SUM", opensn::FieldFunctionInterpolationOperation::OP_SUM)
      .addVariable("OP_AVG", opensn::FieldFunctionInterpolationOperation::OP_AVG)
      .addVariable("OP_MAX", opensn::FieldFunctionInterpolationOperation::OP_MAX)
      .addVariable("OP_SUM_FUNC", opensn::FieldFunctionInterpolationOperation::OP_SUM_FUNC)
      .addVariable("OP_AVG_FUNC", opensn::FieldFunctionInterpolationOperation::OP_AVG_FUNC)
      .addVariable("OP_MAX_FUNC", opensn::FieldFunctionInterpolationOperation::OP_MAX_FUNC);

    luabridge::getGlobalNamespace(L)
      .addVariable("GAUSS_LEGENDRE", opensn::ProductQuadratureType::GAUSS_LEGENDRE)
      .addVariable("GAUSS_CHEBYSHEV", opensn::ProductQuadratureType::GAUSS_CHEBYSHEV)
      .addVariable("GAUSS_LEGENDRE_CHEBYSHEV",
                   opensn::ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV)
      .addVariable("CUSTOM_QUADRATURE", opensn::ProductQuadratureType::CUSTOM_QUADRATURE);
    //

    luabridge::getGlobalNamespace(L)
      .beginClass<LuaLogger>("log")
      .addStaticFunction("Log", &LuaLogger::Log)
      .addStaticFunction("Log0", &LuaLogger::Log0)
      .addStaticFunction("Log0Warning", &LuaLogger::Log0Warning)
      .addStaticFunction("Log0Error", &LuaLogger::Log0Error)
      .addStaticFunction("Log0Verbose0", &LuaLogger::Log0Verbose0)
      .addStaticFunction("Log0Verbose1", &LuaLogger::Log0Verbose1)
      .addStaticFunction("Log0Verbose2", &LuaLogger::Log0Verbose2)
      .addStaticFunction("LogAll", &LuaLogger::LogAll)
      .addStaticFunction("LogAllWarning", &LuaLogger::LogAllWarning)
      .addStaticFunction("LogAllError", &LuaLogger::LogAllError)
      .addStaticFunction("LogAllVerbose0", &LuaLogger::LogAllVerbose0)
      .addStaticFunction("LogAllVerbose1", &LuaLogger::LogAllVerbose1)
      .addStaticFunction("LogAllVerbose2", &LuaLogger::LogAllVerbose2)
      .endClass();

    luabridge::getGlobalNamespace(L).beginNamespace("math").endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginClass<Vector3>("Vector3")
      .addConstructor<void (*)()>()
      .addConstructor<void (*)(double)>()
      .addConstructor<void (*)(double, double)>()
      .addConstructor<void (*)(double, double, double)>()
      .addProperty("x", &Vector3::x)
      .addProperty("y", &Vector3::y)
      .addProperty("z", &Vector3::z)
      .endClass()
      .beginClass<std::shared_ptr<Vector3>>("Vector3Ptr")
      .endClass();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("squad")
      .beginClass<GaussQuadrature>("GaussQuadrature")
      .addProperty("qpoints", &GaussQuadrature::qpoints)
      .addProperty("weights", &GaussQuadrature::weights)
      .endClass()
      .deriveClass<GaussLegendreQuadrature, GaussQuadrature>("GaussLegendreQuadrature")
      .addStaticFunction("Create", &GaussLegendreQuadrature::Create)
      .endClass()
      .beginClass<std::shared_ptr<GaussLegendreQuadrature>>("GaussLegendreQuadraturePtr")
      .endClass()
      .deriveClass<GaussChebyshevQuadrature, GaussQuadrature>("GaussChebyshevQuadrature")
      .addStaticFunction("Create", &GaussChebyshevQuadrature::Create)
      .endClass()
      .beginClass<std::shared_ptr<GaussChebyshevQuadrature>>("GaussChebyshevQuadraturePtr")
      .endClass()
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("aquad")
      .addFunction("Legendre", &Legendre)
      .addFunction("LegendreDerivative", &dLegendredx)
      .addFunction("Ylm", &Ylm)
      .addFunction("CreateProductQuadrature", &AQuadCreateProductQuadrature)
      .addFunction("CreateCylindricalProductQuadrature", &AQuadCreateCylindricalProductQuadrature)
      .addFunction("OptimizeForPolarSymmetry", &AQuadOptimizeForPolarSymmetry)
      .addFunction("CreateSLDFESQuadrature", &AQuadCreateSLDFESQAngularQuadrature)
      .addFunction("LocallyRefineSLDFESQ", &AQuadLocallyRefineSLDFESQ)
      .addFunction("PrintQuadratureToFile", &AQuadPrintQuadratureToFile)
      //
      .beginClass<AngularQuadrature>("AngularQuadrature")
      .endClass()
      .beginClass<std::shared_ptr<AngularQuadrature>>("AngularQuadraturePtr")
      .endClass()
      //
      .deriveClass<ProductQuadrature, AngularQuadrature>("ProductQuadrature")
      .endClass()
      .beginClass<std::shared_ptr<ProductQuadrature>>("ProductQuadraturePtr")
      .endClass()
      //
      .deriveClass<AngularQuadratureProdGL, ProductQuadrature>("AngularQuadratureProdGL")
      .endClass()
      .beginClass<std::shared_ptr<AngularQuadratureProdGL>>("AngularQuadratureProdGLPtr")
      .endClass()
      .deriveClass<SimplifiedLDFESQ::Quadrature, AngularQuadrature>("SLDFESQuadrature")
      .endClass()
      .beginClass<std::shared_ptr<SimplifiedLDFESQ::Quadrature>>("SLDFESQuadraturePtr")
      .endClass()
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginClass<Function>("Function")
      .endClass()
      .beginClass<std::shared_ptr<Function>>("FunctionPtr")
      .endClass()
      //
      .deriveClass<ScalarMaterialFunction, Function>("ScalarMaterialFunction")
      .addFunction("Evaluate", &ScalarMaterialFunction::Evaluate)
      .endClass()
      .beginClass<std::shared_ptr<ScalarMaterialFunction>>("ScalarMaterialFunctionPtr")
      .endClass()
      //
      .deriveClass<LuaScalarMaterialFunction, ScalarMaterialFunction>("LuaScalarMaterialFunction")
      .addStaticFunction("Create", &LuaScalarMaterialFunction::Create)
      .endClass()
      .beginClass<std::shared_ptr<LuaScalarMaterialFunction>>("LuaScalarMaterialFunctionPtr")
      .endClass()
      //
      .deriveClass<ScalarSpatialFunction, Function>("ScalarSpatialFunction")
      .addFunction("Evaluate", &ScalarSpatialFunction::Evaluate)
      .endClass()
      .beginClass<std::shared_ptr<ScalarSpatialFunction>>("ScalarSpatialFunctionPtr")
      .endClass()
      //
      .deriveClass<LuaScalarSpatialFunction, ScalarSpatialFunction>("LuaScalarSpatialFunction")
      .addStaticFunction("Create", &LuaScalarSpatialFunction::Create)
      .endClass()
      .beginClass<std::shared_ptr<LuaScalarSpatialFunction>>("LuaScalarSpatialFunctionPtr")
      .endClass()
      //
      .deriveClass<ScalarSpatialMaterialFunction, Function>("ScalarSpatialMaterialFunction")
      .addFunction("Evaluate", &ScalarSpatialMaterialFunction::Evaluate)
      .endClass()
      .beginClass<std::shared_ptr<ScalarSpatialMaterialFunction>>(
        "ScalarSpatialMaterialFunctionPtr")
      .endClass()
      //
      .deriveClass<LuaScalarSpatialMaterialFunction, ScalarSpatialMaterialFunction>(
        "LuaScalarSpatialMaterialFunction")
      .addStaticFunction("Create", &LuaScalarSpatialMaterialFunction::Create)
      .endClass()
      .beginClass<std::shared_ptr<LuaScalarSpatialMaterialFunction>>(
        "LuaScalarSpatialMaterialFunctionPtr")
      .endClass()
      //
      .deriveClass<VectorSpatialFunction, Function>("VectorSpatialFunction")
      .addFunction("Evaluate", &VectorSpatialFunction::Evaluate)
      .endClass()
      .beginClass<std::shared_ptr<VectorSpatialFunction>>("VectorSpatialFunctionPtr")
      .endClass()
      //
      .deriveClass<LuaVectorSpatialFunction, VectorSpatialFunction>("LuaVectorSpatialFunction")
      .addStaticFunction("Create", &LuaVectorSpatialFunction::Create)
      .endClass()
      .beginClass<std::shared_ptr<LuaVectorSpatialFunction>>("LuaVectorSpatialFunctionPtr")
      .endClass()
      //
      .deriveClass<VectorSpatialMaterialFunction, Function>("VectorSpatialMaterialFunction")
      .addFunction("Evaluate", &VectorSpatialMaterialFunction::Evaluate)
      .endClass()
      .beginClass<std::shared_ptr<VectorSpatialMaterialFunction>>(
        "VectorSpatialMaterialFunctionPtr")
      .endClass()
      //
      .deriveClass<LuaVectorSpatialMaterialFunction, VectorSpatialMaterialFunction>(
        "LuaVectorSpatialMaterialFunction")
      .addStaticFunction("Create", &LuaVectorSpatialMaterialFunction::Create)
      .endClass()
      .beginClass<std::shared_ptr<LuaVectorSpatialMaterialFunction>>(
        "LuaVectorSpatialMaterialFunctionPtr")
      .endClass();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("logvol")
      .beginClass<LogicalVolume>("LogicalVolume")
      .addFunction("Inside", &LogicalVolume::Inside)
      .endClass()
      .beginClass<std::shared_ptr<LogicalVolume>>("LogicalVolumePtr")
      .endClass()
      //
      .deriveClass<BooleanLogicalVolume, LogicalVolume>("BooleanLogicalVolume")
      .addStaticFunction("Create", &BooleanLogicalVolume::Create)
      .endClass()
      .beginClass<std::shared_ptr<BooleanLogicalVolume>>("BooleanLogicalVolumePtr")
      .endClass()
      //
      .deriveClass<RCCLogicalVolume, LogicalVolume>("RCCLogicalVolume")
      .addStaticFunction("Create", &RCCLogicalVolume::Create)
      .endClass()
      .beginClass<std::shared_ptr<RCCLogicalVolume>>("RCCLogicalVolumePtr")
      .endClass()
      //
      .deriveClass<RPPLogicalVolume, LogicalVolume>("RPPLogicalVolume")
      .addStaticFunction("Create", &RPPLogicalVolume::Create)
      .endClass()
      .beginClass<std::shared_ptr<RPPLogicalVolume>>("RPPLogicalVolumePtr")
      .endClass()
      //
      .deriveClass<SphereLogicalVolume, LogicalVolume>("SphereLogicalVolume")
      .addStaticFunction("Create", &SphereLogicalVolume::Create)
      .endClass()
      .beginClass<std::shared_ptr<SphereLogicalVolume>>("SphereLogicalVolumePtr")
      .endClass()
      //
      .deriveClass<SurfaceMeshLogicalVolume, LogicalVolume>("SurfaceMeshLogicalVolume")
      .addStaticFunction("Create", &SurfaceMeshLogicalVolume::Create)
      .endClass()
      .beginClass<std::shared_ptr<SurfaceMeshLogicalVolume>>("SurfaceMeshLogicalVolumePtr")
      .endClass()
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("mesh")
      .addFunction("SetUniformMaterialID", &MeshSetUniformMaterialID)
      .addFunction("SetMaterialIDFromLogicalVolume", &MeshSetMaterialIDFromLogicalVolume)
      .addFunction("SetBoundaryIDFromLogicalVolume", &MeshSetBoundaryIDFromLogicalVolume)
      .addFunction("SetupOrthogonalBoundaries", &MeshSetupOrthogonalBoundaries)
      .addFunction("SetMaterialIDFromFunction", &MeshSetMaterialIDFromFunction)
      .addFunction("SetBoundaryIDFromFunction", &MeshSetBoundaryIDFromFunction)
      .addFunction("ExportToPVTU", &MeshExportToPVTU)
      .addFunction("ComputeVolumePerMaterialID", &MeshComputeVolumePerMaterialID)
      .beginClass<MeshContinuum>("MeshContinuum")
      .addFunction("Dimension", &MeshContinuum::GetDimension)
      .addFunction("SetDimension", &MeshContinuum::SetDimension)
      .endClass()
      .beginClass<std::shared_ptr<MeshContinuum>>("MeshContinuumPtr")
      .endClass()
      //
      .beginClass<MeshGenerator>("MeshGenerator")
      .addStaticFunction("Create", &MeshGenerator::Create)
      .addFunction("Execute", &MeshGenerator::Execute)
      .endClass()
      .beginClass<std::shared_ptr<MeshGenerator>>("MeshGeneratorPtr")
      .endClass()
      //
      .deriveClass<ExtruderMeshGenerator, MeshGenerator>("ExtruderMeshGenerator")
      .addStaticFunction("Create", &ExtruderMeshGenerator::Create)
      .endClass()
      .beginClass<std::shared_ptr<ExtruderMeshGenerator>>("ExtruderMeshGeneratorPtr")
      .endClass()
      //
      .deriveClass<OrthogonalMeshGenerator, MeshGenerator>("OrthogonalMeshGenerator")
      .addStaticFunction("Create", &OrthogonalMeshGenerator::Create)
      .endClass()
      .beginClass<std::shared_ptr<OrthogonalMeshGenerator>>("OrthogonalMeshGeneratorPtr")
      .endClass()
      //
      .deriveClass<FromFileMeshGenerator, MeshGenerator>("FromFileMeshGenerator")
      .addStaticFunction("Create", &FromFileMeshGenerator::Create)
      .endClass()
      .beginClass<std::shared_ptr<FromFileMeshGenerator>>("FromFileMeshGeneratorPtr")
      .endClass()
      //
      .deriveClass<SplitFileMeshGenerator, MeshGenerator>("SplitFileMeshGenerator")
      .addStaticFunction("Create", &SplitFileMeshGenerator::Create)
      .endClass()
      .beginClass<std::shared_ptr<SplitFileMeshGenerator>>("SplitFileMeshGeneratorPtr")
      .endClass()
      //
      .deriveClass<DistributedMeshGenerator, MeshGenerator>("DistributedMeshGenerator")
      .addStaticFunction("Create", &DistributedMeshGenerator::Create)
      .endClass()
      .beginClass<std::shared_ptr<DistributedMeshGenerator>>("DistributedMeshGeneratorPtr")
      .endClass()
      //
      .beginClass<GraphPartitioner>("GraphPartitioner")
      .endClass()
      //
      .deriveClass<KBAGraphPartitioner, GraphPartitioner>("KBAGraphPartitioner")
      .addStaticFunction("Create", &KBAGraphPartitioner::Create)
      .endClass()
      .beginClass<std::shared_ptr<KBAGraphPartitioner>>("KBAGraphPartitionerPtr")
      .endClass()
      //
      .deriveClass<LinearGraphPartitioner, GraphPartitioner>("LinearGraphPartitioner")
      .addStaticFunction("Create", &LinearGraphPartitioner::Create)
      .endClass()
      .beginClass<std::shared_ptr<LinearGraphPartitioner>>("LinearGraphPartitionerPtr")
      .endClass()
      //
      .deriveClass<PETScGraphPartitioner, GraphPartitioner>("PETScGraphPartitioner")
      .addStaticFunction("Create", &PETScGraphPartitioner::Create)
      .endClass()
      .beginClass<std::shared_ptr<PETScGraphPartitioner>>("PETScGraphPartitionerPtr")
      .endClass()
      //
      .beginClass<SurfaceMesh>("SurfaceMesh")
      .addStaticFunction("Create", &SurfaceMesh::Create)
      .addFunction("ImportFromOBJFile", &SurfaceMesh::ImportFromOBJFile)
      .addFunction("ImportFromTriangleFiles", &SurfaceMesh::ImportFromTriangleFiles)
      .addFunction("ImportFromMshFiles", &SurfaceMesh::ImportFromMshFiles)
      .endClass()
      .beginClass<std::shared_ptr<SurfaceMesh>>("SurfaceMeshPtr")
      .endClass()
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("fieldfunc")
      .beginClass<FieldFunction>("FieldFunction")
      .endClass()
      .beginClass<std::shared_ptr<FieldFunction>>("FieldFunctionPtr")
      .endClass()
      //
      .deriveClass<FieldFunctionGridBased, FieldFunction>("FieldFunctionGridBased")
      .endClass()
      .beginClass<std::shared_ptr<FieldFunctionGridBased>>("FieldFunctionGridBasedPtr")
      .endClass()
      //
      .beginClass<FieldFunctionInterpolation>("FieldFunctionInterpolation")
      .addFunction("AddFieldFunction", &FieldFunctionInterpolation::AddFieldFunction)
      .addFunction("Initialize", &FieldFunctionInterpolation::Initialize)
      .addFunction("Execute", &FieldFunctionInterpolation::Execute)
      .addFunction("ExportToCSV", &FieldFunctionInterpolation::ExportToCSV)
      .addFunction("ExportToPython", &FieldFunctionInterpolation::ExportToPython)
      .endClass()
      .beginClass<std::shared_ptr<FieldFunctionInterpolation>>("FieldFunctionInterpolationPtr")
      .endClass()
      //
      .deriveClass<FieldFunctionInterpolationPoint, FieldFunctionInterpolation>(
        "FieldFunctionInterpolationPoint")
      .addFunction("GetPointValue", &FieldFunctionInterpolationPoint::GetPointValue)
      .addStaticFunction("Create", &FieldFunctionInterpolationPoint::Create)
      .endClass()
      .beginClass<std::shared_ptr<FieldFunctionInterpolationPoint>>(
        "FieldFunctionInterpolationPointPtr")
      .endClass()
      //
      .deriveClass<FieldFunctionInterpolationLine, FieldFunctionInterpolation>(
        "FieldFunctionInterpolationLine")
      .addConstructor<void (*)()>()
      .addFunction("SetInitialPoint", &FieldFunctionInterpolationLine::SetInitialPoint)
      .addFunction("SetFinalPoint", &FieldFunctionInterpolationLine::SetFinalPoint)
      .addFunction("SetNumberOfPoints", &FieldFunctionInterpolationLine::SetNumberOfPoints)
      .addStaticFunction("Create", &FieldFunctionInterpolationLine::Create)
      .endClass()
      .beginClass<std::shared_ptr<FieldFunctionInterpolationLine>>(
        "FieldFunctionInterpolationLinePtr")
      .endClass()
      //
      .deriveClass<FieldFunctionInterpolationVolume, FieldFunctionInterpolation>(
        "FieldFunctionInterpolationVolume")
      .addConstructor<void (*)()>()
      .addFunction("SetLogicalVolume", &FieldFunctionInterpolationVolume::SetLogicalVolume)
      .addFunction("SetOperationType", &FieldFunctionInterpolationVolume::SetOperationType)
      .addFunction("SetOperationFunction", &FieldFunctionInterpolationVolume::SetOperationFunction)
      .addFunction("GetValue", &FieldFunctionInterpolationVolume::GetValue)
      .addStaticFunction("Create", &FieldFunctionInterpolationVolume::Create)
      .endClass()
      .beginClass<std::shared_ptr<FieldFunctionInterpolationVolume>>(
        "FieldFunctionInterpolationVolumePtr")
      .endClass()
      //
      .addFunction("ExportToVTK", &ExportFieldFunctionToVTK)
      .addFunction("ExportToVTKMulti", &ExportMultiFieldFunctionToVTK)
      .addFunction("ExportToCSV", &FFInterpolationExportToCSV)
      .addFunction("GetHandleByName", &GetFieldFunctionHandleByName)
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("mat")
      .beginClass<MaterialProperty>("MaterialProperty")
      .endClass()
      .beginClass<std::shared_ptr<MaterialProperty>>("MaterialPropertyPtr")
      .endClass()
      //
      .beginClass<Material>("Material")
      .addFunction("SetTransportXSections", &Material::SetTransportXSections)
      .addFunction("SetIsotropicMGSource", &Material::SetIsotropicMGSource)
      .endClass()
      .beginClass<std::shared_ptr<Material>>("MaterialPtr")
      .endClass()
      //
      .addFunction("AddMaterial", &MatAddMaterial)
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("xs")
      .deriveClass<MultiGroupXS, MaterialProperty>("MultiGroupXS")
      .addFunction("Initialize", luabridge::overload<double, double>(&MultiGroupXS::Initialize))
      .addFunction("SetScalingFactor", &MultiGroupXS::SetScalingFactor)
      .addProperty("num_groups", &MultiGroupXS::GetNumGroups)
      .addProperty("scattering_order", &MultiGroupXS::GetScatteringOrder)
      .addProperty("num_precursors", &MultiGroupXS::GetNumPrecursors)
      .addProperty("is_fissionable", &MultiGroupXS::IsFissionable)
      .addProperty("scaling_factor", &MultiGroupXS::GetScalingFactor)
      .addProperty("sigma_t", &MultiGroupXS::GetSigmaTotal)
      .addProperty("sigma_a", &MultiGroupXS::GetSigmaAbsorption)
      .addProperty("sigma_f", &MultiGroupXS::GetSigmaFission)
      .addProperty("chi", &MultiGroupXS::GetChi)
      .addProperty("nu_sigma_f", &MultiGroupXS::GetNuSigmaF)
      .addProperty("nu_prompt_sigma_f", &MultiGroupXS::GetNuPromptSigmaF)
      .addProperty("nu_delayed_sigma_f", &MultiGroupXS::GetNuDelayedSigmaF)
      .addProperty("inv_velocity", &MultiGroupXS::GetInverseVelocity)
      .endClass()
      .beginClass<std::shared_ptr<MultiGroupXS>>("MultiGroupXSPtr")
      .endClass()
      //
      .deriveClass<IsotropicMultiGroupSource, MaterialProperty>("IsotropicMultiGroupSource")
      .addStaticFunction("FromArray", &IsotropicMultiGroupSource::FromArray)
      .endClass()
      .beginClass<std::shared_ptr<IsotropicMultiGroupSource>>("IsotropicMultiGroupSourcePtr")
      .endClass()
      //
      .addFunction("Create", &XSCreate)
      .addFunction("CreateSimpleOneGroup", &XSCreateSimpleOneGroup)
      .addFunction("LoadFromOpenSn", &XSLoadFromOpenSn)
      .addFunction("LoadFromOpenMC", &XSLoadFromOpenMC)
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("post")
      .beginClass<PostProcessor>("PostProcessor")
      .addFunction("GetValue", &PostProcessor::GetValue)
      .endClass()
      .beginClass<std::shared_ptr<PostProcessor>>("PostProcessorPtr")
      .endClass()
      //
      .deriveClass<SolverInfoPostProcessor, PostProcessor>("SolverInfoPostProcessor")
      .addStaticFunction("Create", &SolverInfoPostProcessor::Create)
      .endClass()
      .beginClass<std::shared_ptr<SolverInfoPostProcessor>>("SolverInfoPostProcessorPtr")
      .endClass()
      //
      .deriveClass<AggregateNodalValuePostProcessor, PostProcessor>(
        "AggregateNodalValuePostProcessor")
      .addStaticFunction("Create", &AggregateNodalValuePostProcessor::Create)
      .endClass()
      .beginClass<std::shared_ptr<AggregateNodalValuePostProcessor>>(
        "AggregateNodalValuePostProcessorPtr")
      .endClass()
      //
      .deriveClass<CellVolumeIntegralPostProcessor, PostProcessor>(
        "CellVolumeIntegralPostProcessor")
      .addStaticFunction("Create", &CellVolumeIntegralPostProcessor::Create)
      .endClass()
      .beginClass<std::shared_ptr<CellVolumeIntegralPostProcessor>>(
        "CellVolumeIntegralPostProcessorPtr")
      .endClass()
      //
      .addFunction("SetPrinterOptions", &PostProcessorPrinterSetOptions)
      .addFunction("Print", &PrintPostProcessors)
      .addFunction("Execute", &ExecutePostProcessors)
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("solver")
      .addFunction("Initialize", &SolverInitialize)
      .addFunction("Execute", &SolverExecute)
      .addFunction("Step", &SolverStep)
      .addFunction("Advance", &SolverAdvance)
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginClass<Solver>("Solver")
      .addFunction("Initialize", &Solver::Initialize)
      .addFunction("Execute", &Solver::Execute)
      .addFunction("Step", &Solver::Step)
      .addFunction("Advance", &Solver::Advance)
      .addFunction("GetFieldFunctions", luabridge::constOverload<>(&Solver::GetFieldFunctions))
      .endClass()
      .beginClass<std::shared_ptr<Solver>>("SolverPtr")
      .endClass();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("diffusion")
      .deriveClass<DiffusionSolverBase, Solver>("DiffusionSolverBase")
      .addFunction("UpdateFieldFunctions", &DiffusionSolverBase::UpdateFieldFunctions)
      .endClass()
      //
      .deriveClass<CFEMDiffusionSolver, DiffusionSolverBase>("CFEMDiffusionSolver")
      .addStaticFunction("Create", &CFEMDiffusionSolver::Create)
      .addFunction("SetOptions", &CFEMDiffusionSolver::SetOptions)
      .addFunction("SetBoundaryOptions", &CFEMDiffusionSolver::SetBoundaryOptions)
      .addFunction("SetDCoefFunction", &CFEMDiffusionSolver::SetDCoefFunction)
      .addFunction("SetQExtFunction", &CFEMDiffusionSolver::SetQExtFunction)
      .addFunction("SetSigmaAFunction", &CFEMDiffusionSolver::SetSigmaAFunction)
      .endClass()
      .beginClass<std::shared_ptr<CFEMDiffusionSolver>>("CFEMDiffusionSolverPtr")
      .endClass()
      //
      .deriveClass<DFEMDiffusionSolver, DiffusionSolverBase>("DFEMDiffusionSolver")
      .addStaticFunction("Create", &DFEMDiffusionSolver::Create)
      .addFunction("SetOptions", &DFEMDiffusionSolver::SetOptions)
      .addFunction("SetBoundaryOptions", &DFEMDiffusionSolver::SetBoundaryOptions)
      .addFunction("SetDCoefFunction", &DFEMDiffusionSolver::SetDCoefFunction)
      .addFunction("SetQExtFunction", &DFEMDiffusionSolver::SetQExtFunction)
      .addFunction("SetSigmaAFunction", &DFEMDiffusionSolver::SetSigmaAFunction)
      .endClass()
      .beginClass<std::shared_ptr<DFEMDiffusionSolver>>("DFEMDiffusionSolverPtr")
      .endClass()
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("prk")
      .deriveClass<PRKSolver, Solver>("PRKSolver")
      .addStaticFunction("Create", &PRKSolver::Create)
      .addFunction("PopulationPrev", &PRKSolver::GetPopulationPrev)
      .addFunction("PopulationNew", &PRKSolver::GetPopulationNew)
      .addFunction("Period", &PRKSolver::GetPeriod)
      .addFunction("TimePrev", &PRKSolver::GetTimePrev)
      .addFunction("TimeNew", &PRKSolver::GetTimeNew)
      .addFunction("SolutionPrev", &PRKSolver::GetSolutionPrev)
      .addFunction("SolutionNew", &PRKSolver::GetSolutionNew)
      .addFunction("SetRho", &PRKSolver::SetRho)
      .endClass()
      .beginClass<std::shared_ptr<PRKSolver>>("PRKSolverPtr")
      .endClass()
      .endNamespace();

    luabridge::getGlobalNamespace(L)
      .beginNamespace("lbs")
      .deriveClass<LBSSolver, Solver>("LBSSolver")
      .addFunction("SetOptions", &LBSSolver::SetOptions)
      .endClass()
      //
      .deriveClass<DiscreteOrdinatesSolver, LBSSolver>("DiscreteOrdinatesSolver")
      .addStaticFunction("Create", &DiscreteOrdinatesSolver::Create)
      .addFunction("ComputeBalance", &DiscreteOrdinatesSolver::ComputeBalance)
      .endClass()
      .beginClass<std::shared_ptr<DiscreteOrdinatesSolver>>("DiscreteOrdinatesSolverPtr")
      .endClass()
      //
      .deriveClass<DiscreteOrdinatesCurvilinearSolver, DiscreteOrdinatesSolver>(
        "DiscreteOrdinatesCurvilinearSolver")
      .addStaticFunction("Create", &DiscreteOrdinatesCurvilinearSolver::Create)
      .endClass()
      .beginClass<std::shared_ptr<DiscreteOrdinatesCurvilinearSolver>>(
        "DiscreteOrdinatesCurvilinearSolverPtr")
      .endClass()
      //
      .deriveClass<SteadyStateSolver, Solver>("SteadyStateSolver")
      .addStaticFunction("Create", &SteadyStateSolver::Create)
      .endClass()
      .beginClass<std::shared_ptr<SteadyStateSolver>>("SteadyStateSolverPtr")
      .endClass()
      //
      .deriveClass<NonLinearKEigen, Solver>("NonLinearKEigen")
      .addStaticFunction("Create", &NonLinearKEigen::Create)
      .endClass()
      .beginClass<std::shared_ptr<NonLinearKEigen>>("NonLinearKEigenPtr")
      .endClass()
      //
      .deriveClass<PowerIterationKEigen, Solver>("PowerIterationKEigen")
      .addStaticFunction("Create", &PowerIterationKEigen::Create)
      .endClass()
      .beginClass<std::shared_ptr<PowerIterationKEigen>>("PowerIterationKEigenPtr")
      .endClass()
      //
      .deriveClass<PowerIterationKEigenSCDSA, Solver>("PowerIterationKEigenSCDSA")
      .addStaticFunction("Create", &PowerIterationKEigenSCDSA::Create)
      .endClass()
      .beginClass<std::shared_ptr<PowerIterationKEigenSCDSA>>("PowerIterationKEigenSCDSAPtr")
      .endClass()
      //
      .deriveClass<PowerIterationKEigenSMM, Solver>("PowerIterationKEigenSMM")
      .addStaticFunction("Create", &PowerIterationKEigenSMM::Create)
      .endClass()
      .beginClass<std::shared_ptr<PowerIterationKEigenSMM>>("PowerIterationKEigenSMMPtr")
      .endClass()
      //
      .deriveClass<DiffusionDFEMSolver, LBSSolver>("DiffusionDFEMSolver")
      .addStaticFunction("Create", &DiffusionDFEMSolver::Create)
      .endClass()
      .beginClass<std::shared_ptr<DiffusionDFEMSolver>>("DiffusionDFEMSolverPtr")
      .endClass()
      //
      .beginClass<PointSource>("PointSource")
      .addStaticFunction("Create", &PointSource::Create)
      .endClass()
      .beginClass<std::shared_ptr<PointSource>>("PointSourcePtr")
      .endClass()
      //
      .beginClass<VolumetricSource>("VolumetricSource")
      .addStaticFunction("Create", &VolumetricSource::Create)
      .endClass()
      .beginClass<std::shared_ptr<VolumetricSource>>("VolumetricSourcePtr")
      .endClass()
      //
      .beginClass<ResponseEvaluator>("ResponseEvaluator")
      .addStaticFunction("Create", &ResponseEvaluator::Create)
      .addFunction("EvaluateResponse", &ResponseEvaluator::EvaluateResponse)
      .endClass()
      .beginClass<std::shared_ptr<ResponseEvaluator>>("ResponseEvaluatorPtr")
      .endClass()
      //
      .addFunction("GetScalarFieldFunctionList", &LBSGetScalarFieldFunctionList)
      .addFunction("ComputeLeakage", &LBSComputeLeakage)
      .addFunction("WriteFluxMoments", &LBSWriteFluxMoments)
      .addFunction("CreateAndWriteSourceMoments", &LBSCreateAndWriteSourceMoments)
      .addFunction("ReadFluxMomentsAndMakeSourceMoments", &LBSReadFluxMomentsAndMakeSourceMoments)
      .addFunction("ReadSourceMoments", &LBSReadSourceMoments)
      .addFunction("ReadFluxMoments", &LBSReadFluxMoments)
      .addFunction("WriteAngularFluxes", &LBSWriteAngularFluxes)
      .addFunction("ReadAngularFluxes", &LBSReadAngularFluxes)
      .endNamespace();
  });

} // namespace opensnlua
