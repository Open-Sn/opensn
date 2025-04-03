// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/runtime.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/post_processors/post_processor.h"
#include "framework/event_system/event.h"
#include "framework/math/math.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "config.h"
#include "caliper/cali.h"
#include "hdf5.h"
#include <iostream>

namespace opensn
{

// Global variables
Logger& log = Logger::GetInstance();
mpi::Communicator mpi_comm;
bool use_caliper = false;
std::string cali_config("runtime-report(calc.inclusive=true),max_column_width=80");
cali::ConfigManager cali_mgr;
Timer program_timer;
bool suppress_color = false;
std::filesystem::path input_path;

std::vector<std::shared_ptr<FieldFunctionInterpolation>> field_func_interpolation_stack;
std::vector<std::shared_ptr<MultiGroupXS>> multigroup_xs_stack;
std::vector<std::shared_ptr<FieldFunction>> field_function_stack;
std::vector<std::shared_ptr<Object>> object_stack;
std::vector<std::shared_ptr<SpatialDiscretization>> sdm_stack;
std::vector<std::shared_ptr<PostProcessor>> postprocessor_stack;

int
Initialize()
{
  if (use_caliper)
  {
    cali_mgr.add(cali_config.c_str());
    cali_set_global_string_byname("opensn.version", GetVersionStr().c_str());
    cali_set_global_string_byname("opensn.input", input_path.c_str());
    cali_mgr.start();
  }

  CALI_MARK_BEGIN(opensn::program.c_str());

  SystemWideEventPublisher::GetInstance().PublishEvent(Event("ProgramStart"));

  // Disable internal HDF error reporting
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  return 0;
}

void
Finalize()
{
  SystemWideEventPublisher::GetInstance().PublishEvent(Event("ProgramExecuted"));

  field_func_interpolation_stack.clear();
  multigroup_xs_stack.clear();
  field_function_stack.clear();
  object_stack.clear();
  sdm_stack.clear();
  postprocessor_stack.clear();

  CALI_MARK_END(opensn::program.c_str());
}

void
Exit(int error_code)
{
  mpi_comm.abort(error_code);
}

std::string
GetVersionStr()
{
  return PROJECT_VERSION;
}

} // namespace opensn
