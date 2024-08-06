// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/lbs_iterative_methods.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_group.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/utils/utils.h"
#include "framework/object.h"

namespace opensn
{

class DiffusionMIPSolver;
class LBSSolver;

/**Group set functioning as a collection of groups*/
class LBSGroupset : public Object
{
public:
  int id_;
  std::vector<LBSGroup> groups_;
  std::shared_ptr<AngularQuadrature> quadrature_;
  std::shared_ptr<AngleAggregation> angle_agg_;
  UniqueSOGroupings unique_so_groupings_;
  DirIDToSOMap dir_id_to_so_map_;

  int master_num_grp_subsets_;
  int master_num_ang_subsets_;

  std::vector<SubSetInfo> grp_subset_infos_;

  IterativeMethod iterative_method_;
  AngleAggregationType angleagg_method_;
  double residual_tolerance_;
  int max_iterations_;
  int gmres_restart_intvl_;

  bool allow_cycles_;

  bool apply_wgdsa_;
  bool apply_tgdsa_;
  int wgdsa_max_iters_;
  int tgdsa_max_iters_;
  double wgdsa_tol_;
  double tgdsa_tol_;
  bool wgdsa_verbose_;
  bool tgdsa_verbose_;
  std::string wgdsa_string_;
  std::string tgdsa_string_;

  std::shared_ptr<DiffusionMIPSolver> wgdsa_solver_ = nullptr;
  std::shared_ptr<DiffusionMIPSolver> tgdsa_solver_ = nullptr;

  struct TwoGridAccelerationInfo
  {
    std::map<int, TwoGridCollapsedInfo> map_mat_id_2_tginfo;
    EnergyCollapseScheme scheme = EnergyCollapseScheme::JFULL;
  } tg_acceleration_info_;

  UnknownManager psi_uk_man_;

public:
  static InputParameters GetInputParameters();

  /**Input parameters based constructor.*/
  explicit LBSGroupset(const InputParameters& params, int id, const LBSSolver& lbs_solver);

  explicit LBSGroupset(int id);

  LBSGroupset();

  /**Computes the discrete to moment operator.*/
  void BuildDiscMomOperator(unsigned int scattering_order, GeometryType geometry_type);

  /**Computes the moment to discrete operator.*/
  void BuildMomDiscOperator(unsigned int scattering_order, GeometryType geometry_type);

  /**Constructs the groupset subsets.*/
  void BuildSubsets();

  /**Constructs the groupset subsets.*/
  void PrintSweepInfoFile(size_t ev_tag, const std::string& file_name);

private:
  void Init(int id);
};

} // namespace opensn
