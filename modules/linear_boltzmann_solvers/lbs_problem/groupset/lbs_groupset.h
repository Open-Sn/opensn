// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/linear_solver/linear_system_solver.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/utils/utils.h"

namespace opensn
{

class DiffusionMIPSolver;
class LBSProblem;

/**
 * Shared groupset state.
 *
 * This class is not specific to LBSProblem. It is shared infrastructure used across
 * problem implementations, including discrete ordinates, to carry both generic
 * group-iteration settings and groupset-scoped transport data.
 *
 * This class remains the sole groupset representation used by the transport
 * solvers, so method-specific groupset data continues to live here.
 */
class LBSGroupset
{
public:
  explicit LBSGroupset(const InputParameters& params, int id, const LBSProblem& lbs_problem);

  explicit LBSGroupset(int id);

  LBSGroupset();

  void PrintSweepInfoFile(size_t ev_tag, const std::string& file_name);

  /// Initialize carrier for copying quadrature data to GPU.
  void InitializeGPUCarriers();

  /// Delete carrier and deallocate memory on GPU.
  void ResetGPUCarriers();

  ~LBSGroupset();

  /// Return number of energy groups in the groupset.
  unsigned int GetNumGroups() const;

  int id;
  /// First energy-group index in the groupset.
  unsigned int first_group;
  /// Last energy-group index in the groupset.
  unsigned int last_group;

  // Sweep/discrete-ordinates data kept in the shared groupset by design.
  // This avoids duplicating generic groupset data.
  std::shared_ptr<AngularQuadrature> quadrature;
  std::shared_ptr<AngleAggregation> angle_agg;
  UniqueSOGroupings unique_so_groupings;
  DirIDToSOMap dir_id_to_so_map;

  int master_num_ang_subsets;

  LinearSystemSolver::IterativeMethod iterative_method;
  AngleAggregationType angleagg_method;
  double residual_tolerance;
  unsigned int max_iterations;
  unsigned int gmres_restart_intvl;

  bool allow_cycles;

  bool apply_wgdsa;
  bool apply_tgdsa;
  unsigned int wgdsa_max_iters;
  unsigned int tgdsa_max_iters;
  double wgdsa_tol;
  double tgdsa_tol;
  bool wgdsa_verbose;
  bool tgdsa_verbose;
  std::string wgdsa_string;
  std::string tgdsa_string;

  void* quad_carrier = nullptr;

  std::shared_ptr<DiffusionMIPSolver> wgdsa_solver = nullptr;
  std::shared_ptr<DiffusionMIPSolver> tgdsa_solver = nullptr;

  struct TwoGridAccelerationInfo
  {
    std::map<unsigned int, TwoGridCollapsedInfo> map_mat_id_2_tginfo;
    EnergyCollapseScheme scheme = EnergyCollapseScheme::JFULL;
  } tg_acceleration_info_;

  UnknownManager psi_uk_man_;

private:
  void Init(int id);

public:
  static InputParameters GetInputParameters();
};

} // namespace opensn
