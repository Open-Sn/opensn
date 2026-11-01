// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "adapters/precice/cell_power.h"
#include <filesystem>
#include <memory>
#include <string>

namespace opensn::coupling
{

/**
 * One-way preCICE participant that repeatedly writes fixed, cell-integrated power for
 * the duration of a coupling scheme.
 *
 * The participant uses OpenSn's communicator. Construction and WriteSteadyCellPower are collective.
 */
class PreciceAdapter
{
public:
  /**
   * Create the participant. `length_scale` multiplies OpenSn mesh coordinates to give coupling-mesh
   * coordinates in the receiver's length unit, for example 0.01 for an OpenSn mesh in centimeters
   * and a receiver in meters. It does not change transferred power.
   */
  PreciceAdapter(std::filesystem::path configuration_file,
                 std::string participant_name = "OpenSn",
                 std::string mesh_name = "OpenSn-Cells",
                 std::string power_data_name = "TransportCellPower",
                 double length_scale = 1.0);
  ~PreciceAdapter();

  PreciceAdapter(const PreciceAdapter&) = delete;
  PreciceAdapter& operator=(const PreciceAdapter&) = delete;
  PreciceAdapter(PreciceAdapter&&) noexcept;
  PreciceAdapter& operator=(PreciceAdapter&&) noexcept;

  /**
   * Register the cell centroids, then write the same power in every coupling window and block until
   * the coupling has completed. A 3D OpenSn mesh requires a 3D coupling mesh, a 2D mesh may use a
   * 2D or 3D coupling mesh, and a 1D mesh, whose coordinate is z, requires a 3D coupling mesh. May
   * be called once per adapter.
   */
  void WriteSteadyCellPower(const CellPowerData& cell_data);

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace opensn::coupling
