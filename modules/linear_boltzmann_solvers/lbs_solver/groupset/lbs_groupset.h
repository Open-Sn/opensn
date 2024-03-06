#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/lbs_iterative_methods.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_group.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include "framework/math/quadratures/legendre_poly/legendrepoly.h"
#include "framework/math/quadratures/angular_quadrature_base.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/physics/physics_namespace.h"
#include "framework/utils/utils.h"
#include "framework/object.h"

namespace opensn
{
namespace lbs
{

class DiffusionMIPSolver;
class LBSSolver;

/**Group set functioning as a collection of groups*/
class LBSGroupset : public Object
{
public:
  int id_;
  std::vector<LBSGroup> groups_;
  std::shared_ptr<AngularQuadrature> quadrature_ = nullptr;
  std::shared_ptr<AngleAggregation> angle_agg_;
  UniqueSOGroupings unique_so_groupings_;
  DirIDToSOMap dir_id_to_so_map_;

  int master_num_grp_subsets_ = 1;
  int master_num_ang_subsets_ = 1;

  std::vector<SubSetInfo> grp_subset_infos_;

  IterativeMethod iterative_method_ = IterativeMethod::CLASSICRICHARDSON;
  AngleAggregationType angleagg_method_ = AngleAggregationType::POLAR;
  double residual_tolerance_ = 1.0e-6;
  int max_iterations_ = 200;
  int gmres_restart_intvl_ = 30;

  bool allow_cycles_ = false;
  bool log_sweep_events_ = false;

  bool apply_wgdsa_ = false;
  bool apply_tgdsa_ = false;
  int wgdsa_max_iters_ = 30;
  int tgdsa_max_iters_ = 30;
  double wgdsa_tol_ = 1.0e-4;
  double tgdsa_tol_ = 1.0e-4;
  bool wgdsa_verbose_ = false;
  bool tgdsa_verbose_ = false;
  std::string wgdsa_string_;
  std::string tgdsa_string_;

  std::shared_ptr<DiffusionMIPSolver> wgdsa_solver_;
  std::shared_ptr<DiffusionMIPSolver> tgdsa_solver_;

  struct TwoGridAccelerationInfo
  {
    std::map<int, TwoGridCollapsedInfo> map_mat_id_2_tginfo;
    EnergyCollapseScheme scheme = EnergyCollapseScheme::JFULL;
  } tg_acceleration_info_;

  UnknownManager psi_uk_man_;

  // lbs_groupset.cc
  static InputParameters GetInputParameters();
  /**Input parameters based constructor.*/
  explicit LBSGroupset(const InputParameters& params, int id, const LBSSolver& lbs_solver);
  LBSGroupset() : LBSGroupset(-1){};
  explicit LBSGroupset(int id) : id_(id) {}

  /**Computes the discrete to moment operator.*/
  void BuildDiscMomOperator(unsigned int scattering_order, GeometryType geometry_type);
  /**Computes the moment to discrete operator.*/
  void BuildMomDiscOperator(unsigned int scattering_order, GeometryType geometry_type);
  /**Constructs the groupset subsets.*/
  void BuildSubsets();

public:
  /**Constructs the groupset subsets.*/
  void PrintSweepInfoFile(size_t ev_tag, const std::string& file_name);
};

} // namespace lbs
} // namespace opensn
