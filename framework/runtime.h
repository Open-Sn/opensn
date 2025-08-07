// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpicpp-lite/mpicpp-lite.h"
#include "caliper/cali-manager.h"
#include <utility>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>
#include <filesystem>

namespace mpi = mpicpp_lite;

namespace opensn
{

static const std::string program = "OpenSn";

class FieldFunctionInterpolation;
class Solver;
class MultiGroupXS;
class FieldFunction;
class SpatialDiscretization;
class Timer;
class Logger;

extern mpi::Communicator mpi_comm;
extern Logger& log;
extern Timer program_timer;
extern bool use_caliper;
extern std::string cali_config;
extern cali::ConfigManager cali_mgr;
extern bool suppress_color;
extern std::filesystem::path input_path;

/// Global stack of handlers
extern std::vector<std::shared_ptr<FieldFunctionInterpolation>> field_func_interpolation_stack;
extern std::vector<std::shared_ptr<MultiGroupXS>> multigroup_xs_stack;
extern std::vector<std::shared_ptr<FieldFunction>> field_function_stack;
extern std::vector<std::shared_ptr<SpatialDiscretization>> sdm_stack;

/// Initializes all necessary items
int Initialize();

/// Finalize the run
void Finalize();

/// Gets the version string.
std::string GetVersionStr();

} // namespace opensn
