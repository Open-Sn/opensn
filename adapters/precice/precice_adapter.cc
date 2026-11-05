// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "adapters/precice/precice_adapter.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <precice/precice.hpp>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace opensn::coupling
{
namespace
{
double
ValidatedLengthScale(const double length_scale)
{
  if (not std::isfinite(length_scale) or length_scale <= 0.0)
    throw std::invalid_argument("length_scale must be finite and positive");
  return length_scale;
}
} // namespace

class PreciceAdapter::Impl
{
public:
  Impl(std::filesystem::path configuration_file,
       std::string participant_name,
       std::string mesh_name,
       std::string power_data_name,
       const double length_scale)
    : mesh_name_(std::move(mesh_name)),
      power_data_name_(std::move(power_data_name)),
      length_scale_(ValidatedLengthScale(length_scale)),
      communicator_(mpi_comm),
      participant_(participant_name,
                   configuration_file.string(),
                   mpi_comm.rank(),
                   mpi_comm.size(),
                   static_cast<void*>(&communicator_))
  {
    if (participant_.getDataDimensions(mesh_name_, power_data_name_) != 1)
      throw std::invalid_argument("OpenSn preCICE cell-power data must be scalar");
  }

  void WriteSteadyCellPower(const CellPowerData& cell_data)
  {
    if (used_)
      throw std::logic_error(
        "PreciceAdapter::WriteSteadyCellPower may be called only once per adapter");

    const auto& cells = cell_data.GetCells();
    const double global_power = GetValidatedGlobalPower(cell_data);
    const unsigned int dimension = cell_data.GetDimension();
    const int mesh_dimension = participant_.getMeshDimensions(mesh_name_);
    if (dimension != 2 and mesh_dimension != 3)
      throw std::invalid_argument("a " + std::to_string(dimension) +
                                  "D OpenSn mesh requires a 3D preCICE coupling mesh");
    log.Log() << "OpenSn preCICE sending cell-integrated power = " << global_power << " W";

    std::vector<double> coordinates;
    std::vector<double> powers;
    coordinates.reserve(static_cast<std::size_t>(mesh_dimension) * cells.size());
    powers.reserve(cells.size());
    for (const auto& cell : cells)
    {
      coordinates.push_back(length_scale_ * cell.centroid.x);
      coordinates.push_back(length_scale_ * cell.centroid.y);
      if (mesh_dimension == 3)
        coordinates.push_back(length_scale_ * cell.centroid.z);
      powers.push_back(cell.power);
    }

    used_ = true;
    std::vector<precice::VertexID> vertex_ids(cells.size());
    participant_.setMeshVertices(mesh_name_, coordinates, vertex_ids);

    if (participant_.requiresInitialData())
      participant_.writeData(mesh_name_, power_data_name_, vertex_ids, powers);

    participant_.initialize();
    // The power is fixed, so implicit-coupling checkpoints hold no state; the queries are still
    // required by preCICE.
    while (participant_.isCouplingOngoing())
    {
      participant_.requiresWritingCheckpoint();
      const double time_step = participant_.getMaxTimeStepSize();
      if (not std::isfinite(time_step) or time_step <= 0.0)
        throw std::runtime_error("preCICE returned a nonpositive or non-finite time-step size");
      participant_.writeData(mesh_name_, power_data_name_, vertex_ids, powers);
      participant_.advance(time_step);
      participant_.requiresReadingCheckpoint();
    }
    participant_.finalize();
  }

private:
  std::string mesh_name_;
  std::string power_data_name_;
  double length_scale_;
  // preCICE keeps a pointer to the communicator, so it must outlive participant_.
  MPI_Comm communicator_;
  precice::Participant participant_;
  bool used_ = false;
};

PreciceAdapter::PreciceAdapter(std::filesystem::path configuration_file,
                               std::string participant_name,
                               std::string mesh_name,
                               std::string power_data_name,
                               const double length_scale)
  : impl_(std::make_unique<Impl>(std::move(configuration_file),
                                 std::move(participant_name),
                                 std::move(mesh_name),
                                 std::move(power_data_name),
                                 length_scale))
{
}

PreciceAdapter::~PreciceAdapter() = default;
PreciceAdapter::PreciceAdapter(PreciceAdapter&&) noexcept = default;
PreciceAdapter& PreciceAdapter::operator=(PreciceAdapter&&) noexcept = default;

void
PreciceAdapter::WriteSteadyCellPower(const CellPowerData& cell_data)
{
  impl_->WriteSteadyCellPower(cell_data);
}

} // namespace opensn::coupling
