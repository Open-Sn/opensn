// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
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
 * configurations of interest in the input. With this, the user can now iterate over the source
 * configurations in the input and convolve them against all available adjoint solutions in the
 * buffer. For example, \code buffers = {
 *        {
 *            buffer1,
 *            buffer2,
 *            ...,
 *            bufferN
 *        }
 *    }
 *
 *    sources = {
 *        {
 *            source1,
 *            source2,
 *            ...,
 *            sourceM
 *        }
 *    }
 *
 *    evaluator = lbs.ResponseEvaluator.Create({ lbs_solver_handle = phys })
 *    evaluator:AddResponseBuffers(buffers)
 *
 *    responses = {}
 *    for i = 1, M do
 *        evaluator:ClearForwardSources()
 *        evaluator:AddResponseSources(sources[i])
 *
 *        responses[i] = {}
 *        for j = 1, N do
 *            responses[i][buffer_name[j]] =
 *                evaluator:EvaluateResponse(buffer_name[j])
 *        end
 *    end
 * \endcode
 */
class ResponseEvaluator
{
private:
  using FluxMomentBuffer = std::vector<double>;
  using AngularFluxBuffer = std::vector<std::vector<double>>;
  using SurfaceAngularFluxBuffer = std::vector<LBSSolverIO::SurfaceAngularFlux>;
  struct AdjointBuffer 
  {
    FluxMomentBuffer flux_moments;
    AngularFluxBuffer angular_fluxes;
    SurfaceAngularFluxBuffer surface_angular_fluxes;
  };
  using MaterialSources = std::map<int, std::vector<double>>;
  using PointSources = std::vector<std::shared_ptr<PointSource>>;
  using VolumetricSources = std::vector<std::shared_ptr<VolumetricSource>>;
  using BoundarySources = std::map<uint64_t, BoundaryPreference>;

public:
  explicit ResponseEvaluator(const InputParameters& params);

  static InputParameters GetOptionsBlock();
  void SetOptions(const InputParameters& params);

  static InputParameters GetBufferOptionsBlock();
  void SetBufferOptions(const InputParameters& input);

  static InputParameters GetSourceOptionsBlock();
  void SetSourceOptions(const InputParameters& input);

  static InputParameters GetMaterialSourceOptionsBlock();
  void SetMaterialSourceOptions(const InputParameters& params);

  void SetBoundarySourceOptions(const InputParameters& params);

  /// Clear the existing forward sources from the response evaluator.
  void ClearForwardSources();

  void AddResponseBuffers(const InputParameters& params);

  void AddResponseSources(const InputParameters& params);

  /**
   * Evaluate a response using the specified adjoint buffer with the currently defined sources in
   * the solver.
   */
  double EvaluateResponse(const std::string& buffer_name) const;

  /**
   * Evaluate a response using the specified forward buffer and 
   * adjoint buffer towards a specified surface id.
   */
  double EvaluateSurfaceResponse(const std::string& fwd_buffer,
                                 const std::string& adj_buffer) const;
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

private:
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
};

} // namespace opensn
