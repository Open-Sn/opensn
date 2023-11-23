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

std::vector<MeshHandlerPtr> Chi::meshhandler_stack;
int Chi::current_mesh_handler = -1;

std::vector<SurfaceMeshPtr> Chi::surface_mesh_stack;
std::vector<FFInterpPtr> Chi::field_func_interpolation_stack;
std::vector<UnpartMeshPtr> Chi::unpartitionedmesh_stack;

std::vector<MaterialPtr> Chi::material_stack;
std::vector<MultiGroupXSPtr> Chi::multigroup_xs_stack;
std::vector<FieldFunctionPtr> Chi::field_function_stack;

std::vector<AngularQuadraturePtr> Chi::angular_quadrature_stack;

std::vector<ChiObjectPtr> Chi::object_stack;
std::vector<SpatialDiscretizationPtr> Chi::sdm_stack;
std::vector<PostProcessorPtr> Chi::postprocessor_stack;
std::vector<FunctionPtr> Chi::function_stack;

bool Chi::suppress_color_ = false;

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
  Chi::meshhandler_stack.clear();

  Chi::surface_mesh_stack.clear();
  Chi::object_stack.clear();
  Chi::field_func_interpolation_stack.clear();
  Chi::unpartitionedmesh_stack.clear();

  Chi::object_stack.clear();
  Chi::material_stack.clear();
  Chi::multigroup_xs_stack.clear();
  Chi::function_stack.clear();
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
