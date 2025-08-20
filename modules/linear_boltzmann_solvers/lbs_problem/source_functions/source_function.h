// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include <memory>
#include <utility>

namespace opensn
{

class LBSProblem;
class LBSGroupset;

/**
 * Implements a customizable source function using virtual methods.
 * This base class will function well for steady simulations and kEigenvalue
 * simulations. It needs some customization for adjoint and transient.
 */
class SourceFunction
{
protected:
  const LBSProblem& lbs_problem_;

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
  /// Constructor.
  explicit SourceFunction(const LBSProblem& lbs_problem);
  virtual ~SourceFunction() = default;

  /**
   * Sets the source moments for the groups in the current group set.
   *
   * \param groupset The groupset the under consideration.
   * \param q A vector to contribute the source to.
   * \param phi The primary STL vector to operate off.
   * \param source_flags Flags for adding specific terms into the
   *        destination vector. Available flags are for applying
   *        the material source, across/within-group scattering,
   *        and across/within-groups fission.
   */
  virtual void operator()(const LBSGroupset& groupset,
                          std::vector<double>& q,
                          const std::vector<double>& phi,
                          SourceFlags source_flags);

  virtual double FixedSourceMoments() const;

  using PrecursorList = std::vector<MultiGroupXS::Precursor>;
  /// Adds delayed particle precursor sources.
  virtual double AddDelayedFission(const PrecursorList& precursors,
                                   const double& rho,
                                   const std::vector<double>& nu_delayed_sigma_f,
                                   const double* phi) const;

  virtual void AddAdditionalSources(const LBSGroupset& groupset,
                                    std::vector<double>& q,
                                    const std::vector<double>& phi,
                                    SourceFlags source_flags)
  {
    AddPointSources(groupset, q, phi, source_flags);
    AddVolumetricSources(groupset, q, phi, source_flags);
  }

  /// Adds point sources to the source moments.
  void AddPointSources(const LBSGroupset& groupset,
                       std::vector<double>& q,
                       const std::vector<double>& phi,
                       SourceFlags source_flags);

  /// Adds volumetric sources to the source moments.
  void AddVolumetricSources(const LBSGroupset& groupset,
                            std::vector<double>& q,
                            const std::vector<double>& phi,
                            SourceFlags source_flags);
};

} // namespace opensn
