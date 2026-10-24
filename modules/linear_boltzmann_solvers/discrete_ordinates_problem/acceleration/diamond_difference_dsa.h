// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/unknown_manager/unknown_manager.h"
#include <mpi.h>
#include <memory>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;
class PieceWiseLinearContinuous;

/**
 * Multigroup DSA consistent with diamond-difference slab transport.
 *
 * The correction of each group is solved on the cell edges in increasing group order with the
 * diamond-difference mass matrix.
 */
class DiamondDifferenceDSA
{
public:
  explicit DiamondDifferenceDSA(DiscreteOrdinatesProblem& problem);
  ~DiamondDifferenceDSA();

  DiamondDifferenceDSA(const DiamondDifferenceDSA&) = delete;
  DiamondDifferenceDSA& operator=(const DiamondDifferenceDSA&) = delete;

  /// Add the scalar-flux correction after a diamond-difference sweep.
  void Apply();

private:
  struct GroupSystem;
  DiscreteOrdinatesProblem& problem_;
  std::shared_ptr<PieceWiseLinearContinuous> edges_;
  UnknownManager unknowns_;
  MPI_Comm communicator_ = MPI_COMM_NULL;
  std::vector<std::unique_ptr<GroupSystem>> systems_;
};

} // namespace opensn
