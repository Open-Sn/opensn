// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/parameters/input_parameters.h"
#include <memory>

namespace opensn
{

/**
 * A class used for evaluating responses by folding sources against adjoint solutions.
 *
 * The workflow for this utility is constructed to minimize the file reading necessary for
 * evaluations. To begin, one should add all adjoint solutions that are desired for response
 * computations into the buffer. Then, one should define the different forward source
 * configurations of interest. With this, the user can now iterate over the source
 * configurations and convolve them against all available adjoint solutions in the buffer. For
 * example, \code{.py}
 *    evaluator = ResponseEvaluator(problem=phys)
 *    evaluator.SetOptions(buffers=[{"name": "buff1", "file_prefixes": {...}}, ...])
 *
 *    responses = {}
 *    for source_config in source_configs:
 *        evaluator.ClearForwardSources()
 *        evaluator.SetOptions(sources=source_config)
 *
 *        responses[source_config] = {
 *            buffer_name: evaluator.EvaluateResponse(buffer_name)
 *            for buffer_name in buffer_names
 *        }
 * \endcode
 */
class ResponseEvaluator
{
private:
  using FluxMomentBuffer = std::vector<double>;
  using AngularFluxBuffer = std::vector<std::vector<double>>;
  using AdjointBuffer = std::pair<FluxMomentBuffer, AngularFluxBuffer>;

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

  void SetBufferOptions(const InputParameters& input);

  void SetSourceOptions(const InputParameters& input);

  void SetMaterialSourceOptions(const InputParameters& params);

  void SetBoundarySourceOptions(const InputParameters& params);

  /// Clear the existing forward sources from the response evaluator.
  void ClearForwardSources();

  void AddResponseBuffers(const InputParameters& params);

  void AddResponseSources(const InputParameters& params);

  static InputParameters GetBoundarySourceOptionsBlock();

  /**
   * Evaluate a response using the specified adjoint buffer with the currently defined sources in
   * the solver.
   */
  double EvaluateResponse(const std::string& buffer_name) const;

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

  std::map<std::string, AdjointBuffer> adjoint_buffers_;

  MaterialSources material_sources_;
  PointSources point_sources_;
  VolumetricSources volumetric_sources_;
  BoundarySources boundary_sources_;

public:
  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();
  static std::shared_ptr<ResponseEvaluator> Create(const ParameterBlock& params);

  /// Return the schema for the options block, covering adjoint buffers and forward sources.
  static InputParameters GetOptionsBlock();
  static InputParameters GetBufferOptionsBlock();
  static InputParameters GetSourceOptionsBlock();
  static InputParameters GetMaterialSourceOptionsBlock();
};

} // namespace opensn
