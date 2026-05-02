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

class Timer;

extern mpi::Communicator mpi_comm;
extern Timer program_timer;
extern bool use_caliper;
extern std::string cali_config;
extern cali::ConfigManager cali_mgr;
extern std::filesystem::path input_path;

/// Initializes all necessary items
int Initialize();

/// Finalize the run
void Finalize();

/// Gets the version string.
std::string GetVersionStr();

} // namespace opensn
