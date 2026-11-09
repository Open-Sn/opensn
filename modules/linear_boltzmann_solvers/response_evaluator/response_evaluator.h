// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/response_evaluator/surface_response_data.h"
#include "framework/parameters/input_parameters.h"
#include <memory>

namespace opensn
{

/**
 * A class used for evaluating responses by folding sources against adjoint solutions.
 *
 * The workflow for this utility is constructed to minimize the file reading necessary for
 * evaluations. To begin, one should add all adjoint solutions that are desired for response
 * computations as named data sets. Then, one should define the different forward source
 * configurations of interest. With this, the user can now iterate over the source
 * configurations and convolve them against all available adjoint data sets. For
 * example, \code{.py}
 *    evaluator = ResponseEvaluator(problem=phys)
 *    evaluator.SetOptions(data=[{"name": "detector", "file_prefixes": {...}}, ...])
 *
 *    responses = {}
 *    for config_index, source_config in enumerate(source_configs):
 *        evaluator.ClearForwardSources()
 *        evaluator.SetOptions(sources=source_config)
 *
 *        responses[config_index] = {
 *            data_name: evaluator.EvaluateResponse(data_name)
 *            for data_name in data_names
 *        }
 * \endcode
 */
class ResponseEvaluator
{
private:
  using FluxMomentData = std::vector<double>;
  using AngularFluxData = std::vector<std::vector<double>>;
  /// Surface angular fluxes distributed by face over the ranks; see SurfaceShard.
  using SurfaceAngularFluxData = DistributedSurfaceData;
  struct FluxData
  {
    FluxMomentData phi;
    AngularFluxData psi;
    SurfaceAngularFluxData surf_psi;
  };

  using MaterialSources = std::map<unsigned int, std::vector<double>>;
  using PointSources = std::vector<std::shared_ptr<PointSource>>;
  using VolumetricSources = std::vector<std::shared_ptr<VolumetricSource>>;
  struct BoundarySource
  {
    LBSBoundaryType type = LBSBoundaryType::VACUUM;
    std::vector<double> isotropic_mg_source;
  };
  using BoundarySources = std::map<uint64_t, BoundarySource>;

public:
  explicit ResponseEvaluator(const InputParameters& params);

  void SetOptions(const InputParameters& params);

  /**
   * Add a named data set. `file_prefixes` may contain `flux_moments` and `angular_fluxes`, the
   * adjoint flux moments and angular fluxes that EvaluateResponse folds against the forward
   * sources, and `surface_angular_fluxes`, the surface data that EvaluateSurfaceResponse pairs.
   * Flux moments and angular fluxes must have been written with the current partitioning. The
   * listed `surfaces` are read from all surface files of the prefix, whatever number of ranks
   * and partitioning wrote them, and distributed by face over the current ranks with
   * DistributeSurfaceAngularFluxes. Each is a boundary name or an interior surface name, which
   * reads both stored sides `<name>_u` and `<name>_d`. Collective.
   */
  void SetDataOptions(const InputParameters& input);

  void SetSourceOptions(const InputParameters& input);

  void SetMaterialSourceOptions(const InputParameters& params);

  void SetBoundarySourceOptions(const InputParameters& params);

  /// Clear the existing forward sources from the response evaluator.
  void ClearForwardSources();

  void AddResponseData(const InputParameters& params);

  void AddResponseSources(const InputParameters& params);

  static InputParameters GetBoundarySourceOptionsBlock();

  /**
   * Evaluate a response by folding the currently defined forward sources against the adjoint flux
   * moments and angular fluxes of the named data set added with SetDataOptions.
   *
   * Collective. Throws `std::invalid_argument` for an unknown data set and `std::logic_error` if
   * the data set lacks the flux moments or angular fluxes that the defined sources require.
   */
  double EvaluateResponse(const std::string& data_name) const;

  /**
   * Evaluate the surface response of forward and adjoint surface angular fluxes added with
   * SetDataOptions.
   *
   * Returns one surface response per selected surface, as defined in ComputeSurfaceResponses: its
   * outgoing and incoming parts and their sum, the response carried across the surface. They are
   * also logged at verbosity 1.
   *
   * \param fwd_data_name Name of the forward data set.
   * \param adj_data_name Name of the adjoint data set.
   * \param surfaces Surfaces to evaluate: boundary names and interior surface names (`<name>` or
   *        `<name>_u` along the `_u` normal, `<name>_d` along the `_d` normal). If empty, every
   *        surface in the forward data is evaluated, in name order.
   *
   * Collective. Throws `std::invalid_argument` for unknown data sets or surfaces, a surface
   * requested more than once, and incompatible data, and `std::logic_error` if the forward and
   * adjoint data do not cover the same faces.
   */
  std::vector<SurfaceResponse>
  EvaluateSurfaceResponse(const std::string& fwd_data_name,
                          const std::string& adj_data_name,
                          const std::vector<std::string>& surfaces = {}) const;

private:
  /**
   * Evaluates a boundary source and returns the angular flux on the boundary.
   *
   * This returns the full angular flux at a particular spatial location, for a particular
   * groupset, at a particular time, for a particular boundary. No boundary normal information
   * is included in the evaluation. The incident fluxes are obtained within the EvaluateResponse
   * routine.
   */
  std::vector<double> EvaluateBoundaryCondition(uint64_t boundary_id,
                                                const Vector3& node,
                                                const LBSGroupset& groupset,
                                                double time = 0.0) const;

  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;

  std::map<std::string, FluxData> flux_data_;

  MaterialSources material_sources_;
  PointSources point_sources_;
  VolumetricSources volumetric_sources_;
  BoundarySources boundary_sources_;

public:
  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();
  static std::shared_ptr<ResponseEvaluator> Create(const ParameterBlock& params);

  /// Return the schema for the options block, covering data sets and forward sources.
  static InputParameters GetOptionsBlock();
  static InputParameters GetDataOptionsBlock();
  static InputParameters GetSourceOptionsBlock();
  static InputParameters GetMaterialSourceOptionsBlock();
};

} // namespace opensn
