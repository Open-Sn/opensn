#ifndef CHITECH_LBS_SOURCE_FUNCTION_H
#define CHITECH_LBS_SOURCE_FUNCTION_H

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_structs.h"

#include "physics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"

#include <memory>
#include <utility>

namespace lbs
{
class LBSSolver;
class LBSGroupset;

//###################################################################
/**Implements a customizable source function using virtual methods.
 * This base class will function well for steady simulations and kEigenvalue
 * simulations. It needs some customization for adjoint and transient.*/
class SourceFunction
{
protected:
  const LBSSolver& lbs_solver_;

  bool apply_fixed_src_ = false;
  bool apply_wgs_scatter_src_ = false;
  bool apply_ags_scatter_src_ = false;
  bool apply_wgs_fission_src_ = false;
  bool apply_ags_fission_src_ = false;
  bool suppress_wg_scatter_src_ = false;

  size_t gs_i_ = 0;
  size_t gs_f_ = 0;
  size_t first_grp_ = 0;
  size_t last_grp_ = 0;

  double cell_volume_ = 0.0;
  size_t g_ = 0;
  const double* fixed_src_moments_ = nullptr;
  std::vector<double> default_zero_src_;

public:
  explicit SourceFunction(const LBSSolver& lbs_solver);
  virtual ~SourceFunction() = default;

  virtual void operator()(LBSGroupset& groupset,
                          std::vector<double>& destination_q,
                          const std::vector<double>& phi,
                          SourceFlags source_flags);

  virtual double AddSourceMoments() const;

  typedef std::vector<chi_physics::MultiGroupXS::Precursor> PrecursorList;
  virtual double AddDelayedFission(const PrecursorList& precursors,
                                   const std::vector<double>& nu_delayed_sigma_f,
                                   const double* phi) const;

  virtual void AddAdditionalSources(LBSGroupset& groupset,
                                    std::vector<double>& destination_q,
                                    const std::vector<double>& phi,
                                    SourceFlags source_flags)
  {
    AddPointSources(groupset, destination_q, phi, source_flags);
  }

  void AddPointSources(LBSGroupset& groupset,
                       std::vector<double>& destination_q,
                       const std::vector<double>& phi,
                       SourceFlags source_flags);
};

} // namespace lbs

#endif // CHITECH_LBS_SOURCE_FUNCTION_H
