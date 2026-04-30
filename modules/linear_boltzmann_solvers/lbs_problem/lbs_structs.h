// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/math.h"
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

enum class PhiSTLOption
{
  PHI_OLD = 1,
  PHI_NEW = 2
};

/// Restart read/write policy and derived per-rank paths.
struct LBSRestartOptions
{
  bool writes_enabled = false;
  bool write_delayed_psi = true;
  std::filesystem::path read_path;
  std::filesystem::path write_path;
  std::chrono::time_point<std::chrono::system_clock> last_write_time;
  std::chrono::seconds write_time_interval = std::chrono::seconds(0);

  bool IsWriteDue() const
  {
    if (write_time_interval <= std::chrono::seconds(0))
      return false;

    const auto elapsed = std::chrono::system_clock::now() - last_write_time;
    return elapsed >= write_time_interval;
  }

  void MarkWriteComplete() { last_write_time = std::chrono::system_clock::now(); }
};

/// Struct for storing LBS options.
struct LBSOptions
{
  int max_mpi_message_size = 32768;

  LBSRestartOptions restart;

  bool use_precursors = true;
  bool use_src_moments = false;
  bool save_angular_flux = false;

  bool adjoint = false;

  bool verbose_inner_iterations = true;
  bool verbose_outer_iterations = true;
  bool ags_pointwise_convergence = false;
  unsigned int max_ags_iterations = 100;
  double ags_tolerance = 1.0e-6;

  double power_default_kappa = 3.20435e-11; // 200MeV to Joule

  std::string field_function_prefix_option = "prefix";
  std::string field_function_prefix; // Default is empty

  LBSOptions() = default;
};

} // namespace opensn
