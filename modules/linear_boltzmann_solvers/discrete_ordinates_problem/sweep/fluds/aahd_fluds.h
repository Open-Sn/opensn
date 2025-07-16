// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/ndarray.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "caribou/caribou.h"
#include <vector>
#include <utility>

namespace crb = caribou;

namespace opensn
{

/**
 * @brief Object storing pointer to group-wise angular flux of a face node and angle on the host.
 * @details If the node is associated with incoming face, the ``outgoing`` member is a null pointer
 * and vice-versa.
 */
struct FloodPtr
{
  const double* incoming = nullptr;
  double* outgoing = nullptr;
};

/// Object storing pointers to angular fluxes for a given face.
struct FloodFace
{
  /**
   * Pointer to fluxes in the FLUDS for each face node (first index) and angle in angleset (second
   * index).
   */
  NDArray<FloodPtr, 2> flux;
};

/// Object storing pointers to angular fluxes for a given cell.
struct FloodCell : public std::vector<FloodFace>
{
};

/// Object storing pointers to angular fluxes for a given SPDS level.
struct FloodLevel
{
  /// Map from cell index in level to the data of the cell faces.
  std::vector<FloodCell> cells;
  /// Precomputed offset index into contiguous memory for incoming flux.
  std::vector<char> precomputed_offset;
};

/** @brief Flatten data structure of AAH_FLUDS for tracking angular flux pointers on a contiguous
 *  memory block on device.
 */
class AAHD_FLUDS
{
public:
  /// @brief Constructor from angle set.
  AAHD_FLUDS(LBSProblem& lbs_problem,
             const LBSGroupset& group_set,
             AngleSet& angle_set,
             bool is_surface_source_active);

  /// Copy data of a given level to GPU.
  void CopyToDevice(int level);

  /// Copy data of a given level back from GPU.
  void CopyFromDevice(int level);

  /// Get the device memory.
  inline char* GetDevicePtr() { return device_buffer_.get(); }

  /// Get strides for each face node.
  inline std::uint32_t GetStride() { return groupset_size_ * angleset_size_; }

protected:
  /**
   * @details Map from level index, index of the cell in the level and face index to the pointer of
   * the corresponding incoming psi in the fluds.
   */
  std::vector<FloodLevel> fluds_map_;
  /// Map from level index to cell details.

  /**
   * @brief Size of required memory for storing incoming flux.
   * @details Minimum size of the contiguous memory required for storing data of the largest level.
   */
  std::size_t size_ = 0;
  /// Number of groups of the associated group set.
  std::size_t groupset_size_ = 0;
  /// Number of angles in the angleset.
  std::size_t angleset_size_ = 0;
  /// Contiguous memory on the host (CPU) for angular flux.
  crb::HostVector<char> host_buffer_;
  /// Contiguous memory on the device (GPU) for angular flux.
  crb::DeviceMemory<char> device_buffer_;
};

} // namespace opensn
