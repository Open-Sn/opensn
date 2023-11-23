/** @file Runtime file*/
#include "framework/runtime.h"
#include "config.h"

#include "framework/math/math.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/physics/physics_namespace.h"

#include "framework/post_processors/post_processor.h"

#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/event_system/event.h"

#include "framework/object_factory.h"

#include "framework/mpi/mpi.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"

#include <iostream>

namespace opensn
{

// Global variables
Logger& log = Logger::GetInstance();
MPI_Info& mpi = MPI_Info::GetInstance();
Timer program_timer;

std::vector<std::shared_ptr<MeshHandler>> meshhandler_stack;
int current_mesh_handler = -1;

std::vector<std::shared_ptr<SurfaceMesh>> surface_mesh_stack;
std::vector<std::shared_ptr<FieldFunctionInterpolation>> field_func_interpolation_stack;
std::vector<std::shared_ptr<UnpartitionedMesh>> unpartitionedmesh_stack;

std::vector<std::shared_ptr<Material>> material_stack;
std::vector<std::shared_ptr<MultiGroupXS>> multigroup_xs_stack;
std::vector<std::shared_ptr<FieldFunction>> field_function_stack;

std::vector<AngularQuadraturePtr> angular_quadrature_stack;

std::vector<std::shared_ptr<Object>> object_stack;
std::vector<SpatialDiscretizationPtr> sdm_stack;
std::vector<PostProcessorPtr> postprocessor_stack;
std::vector<FunctionPtr> function_stack;

bool suppress_color = false;

int
Initialize()
{
  auto& t_main = log.CreateTimingBlock(opensn::name);
  t_main.TimeSectionBegin();
  SystemWideEventPublisher::GetInstance().PublishEvent(Event("ProgramStart"));

  return 0;
}

void
Finalize()
{
  auto& t_main = log.GetTimingBlock(opensn::name);
  t_main.TimeSectionEnd();
  SystemWideEventPublisher::GetInstance().PublishEvent(Event("ProgramExecuted"));
  meshhandler_stack.clear();

  surface_mesh_stack.clear();
  object_stack.clear();
  field_func_interpolation_stack.clear();
  unpartitionedmesh_stack.clear();

  object_stack.clear();
  material_stack.clear();
  multigroup_xs_stack.clear();
  function_stack.clear();
}

void
Exit(int error_code)
{
  MPI_Abort(mpi.comm, error_code);
}

std::string
GetVersionStr()
{
  return PROJECT_VERSION;
}

} // namespace opensn
