#pragma once

#include "framework/ChiObject.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_group.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/lbs_iterativemethods.h"

#include "framework/math/Quadratures/LegendrePoly/legendrepoly.h"
#include "framework/math/Quadratures/angular_quadrature_base.h"
#include "framework/math/UnknownManager/unknown_manager.h"

#include "framework/mesh/SweepUtilities/AngleAggregation/angleaggregation.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

#include "framework/physics/chi_physics_namespace.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/acceleration.h"

#include "framework/utils/chi_utils.h"

namespace lbs::acceleration
{
class DiffusionMIPSolver;
}

namespace lbs
{

class LBSSolver;

// ################################################################### Class def
/**Group set functioning as a collection of groups*/
class LBSGroupset : public ChiObject
{
protected:
  typedef std::shared_ptr<chi_mesh::sweep_management::AngleAggregation> AngleAggPtr;

public:
  int id_;
  std::vector<LBSGroup> groups_;
  std::shared_ptr<chi_math::AngularQuadrature> quadrature_ = nullptr;
  AngleAggPtr angle_agg_;
  UniqueSOGroupings unique_so_groupings_;
  DirIDToSOMap dir_id_to_so_map_;

  int master_num_grp_subsets_ = 1;
  int master_num_ang_subsets_ = 1;

  std::vector<chi::SubSetInfo> grp_subset_infos_;

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

  std::shared_ptr<lbs::acceleration::DiffusionMIPSolver> wgdsa_solver_;
  std::shared_ptr<lbs::acceleration::DiffusionMIPSolver> tgdsa_solver_;

  struct TwoGridAccelerationInfo
  {
    std::map<int, acceleration::TwoGridCollapsedInfo> map_mat_id_2_tginfo;
    acceleration::EnergyCollapseScheme scheme = acceleration::EnergyCollapseScheme::JFULL;
  } tg_acceleration_info_;

  chi_math::UnknownManager psi_uk_man_;

  // lbs_groupset.cc
  static chi::InputParameters GetInputParameters();
  explicit LBSGroupset(const chi::InputParameters& params, int id, const LBSSolver& lbs_solver);
  LBSGroupset() : LBSGroupset(-1){};
  explicit LBSGroupset(int id) : id_(id) {}

  void BuildDiscMomOperator(unsigned int scattering_order, GeometryType geometry_type);
  void BuildMomDiscOperator(unsigned int scattering_order, GeometryType geometry_type);
  void BuildSubsets();

public:
  void PrintSweepInfoFile(size_t ev_tag, const std::string& file_name);
};
} // namespace lbs
