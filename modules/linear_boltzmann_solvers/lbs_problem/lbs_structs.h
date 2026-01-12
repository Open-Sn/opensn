// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/math.h"
#include "framework/math/functions/function.h"
#include <functional>
#include <chrono>
#include <map>
#include <filesystem>
#include <memory>

namespace opensn
{

using BlockID2XSMap = std::map<unsigned int, std::shared_ptr<MultiGroupXS>>;

enum class LBSBoundaryType
{
  VACUUM = 1,     ///< Zero for all angles, space
  ISOTROPIC = 2,  ///< One value for all angles, homogenous in space
  REFLECTING = 3, ///< Reflecting boundary condition about a normal
  ARBITRARY = 4   ///< Complex different for each angle and face node
};

struct BoundaryPreference
{
  LBSBoundaryType type = LBSBoundaryType::VACUUM;
  std::vector<double> isotropic_mg_source;
  std::string source_function;
  std::shared_ptr<AngularFluxFunction> angular_flux_function;
};

enum class PhiSTLOption
{
  PHI_OLD = 1,
  PHI_NEW = 2
};

/// Struct for storing LBS options.
struct LBSOptions
{
  int max_mpi_message_size = 32768;

  bool restart_writes_enabled = false;
  bool write_delayed_psi_to_restart = true;
  std::filesystem::path read_restart_path;
  std::filesystem::path write_restart_path;
  std::chrono::time_point<std::chrono::system_clock> last_restart_write_time;
  std::chrono::seconds write_restart_time_interval = std::chrono::seconds(0);

  bool use_precursors = false;
  bool use_src_moments = false;

  bool save_angular_flux = false;

  bool adjoint = false;

  bool verbose_inner_iterations = true;
  bool verbose_ags_iterations = true;
  bool verbose_outer_iterations = true;
  bool ags_pointwise_convergence = false;
  unsigned int max_ags_iterations = 100;
  double ags_tolerance = 1.0e-6;

  bool power_field_function_on = false;
  double power_default_kappa = 3.20435e-11; // 200MeV to Joule
  double power_normalization = -1.0;

  std::string field_function_prefix_option = "prefix";
  std::string field_function_prefix; // Default is empty

  LBSOptions() = default;
};

} // namespace opensn
