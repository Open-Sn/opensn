// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_group.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/utils/utils.h"
#include "framework/object.h"

namespace opensn
{

class DiffusionMIPSolver;
class LBSProblem;

/// Group set functioning as a collection of groups
class LBSGroupset : public Object
{
public:
  int id;
  std::vector<LBSGroup> groups;
  std::shared_ptr<AngularQuadrature> quadrature;
  std::shared_ptr<AngleAggregation> angle_agg;
  UniqueSOGroupings unique_so_groupings;
  DirIDToSOMap dir_id_to_so_map;

  int master_num_ang_subsets;

  LinearSolver::IterativeMethod iterative_method;
  AngleAggregationType angleagg_method;
  double residual_tolerance;
  int max_iterations;
  int gmres_restart_intvl;

  bool allow_cycles;

  bool apply_wgdsa;
  bool apply_tgdsa;
  int wgdsa_max_iters;
  int tgdsa_max_iters;
  double wgdsa_tol;
  double tgdsa_tol;
  bool wgdsa_verbose;
  bool tgdsa_verbose;
  std::string wgdsa_string;
  std::string tgdsa_string;

  std::shared_ptr<DiffusionMIPSolver> wgdsa_solver = nullptr;
  std::shared_ptr<DiffusionMIPSolver> tgdsa_solver = nullptr;

  struct TwoGridAccelerationInfo
  {
    std::map<int, TwoGridCollapsedInfo> map_mat_id_2_tginfo;
    EnergyCollapseScheme scheme = EnergyCollapseScheme::JFULL;
  } tg_acceleration_info_;

  UnknownManager psi_uk_man_;

public:
  static InputParameters GetInputParameters();

  /// Input parameters based constructor.
  explicit LBSGroupset(const InputParameters& params, int id, const LBSProblem& lbs_problem);

  explicit LBSGroupset(int id);

  LBSGroupset();

  void PrintSweepInfoFile(size_t ev_tag, const std::string& file_name);

  /**
   * Get index of a group within this groupset
   *
   * @param g Group
   * @return Index of the group within this groupset
   */
  std::size_t GetGroupIndex(int g) const;

private:
  void Init(int id);
};

} // namespace opensn
