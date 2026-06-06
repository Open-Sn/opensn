// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/raytrace/raytracer.h"
#include "framework/mesh/cell/cell.h"
#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector.h"
#include "framework/data_types/vector3.h"
#include "hdf5.h"
#include <atomic>
#include <mutex>
#include <set>
#include <string>
#include <thread>

namespace opensn
{

struct UncollidedMatrices
{
  DenseMatrix<double> intV_shapeJ_omega_gradshapeI;
  std::vector<DenseMatrix<double>> intS_omega_n_shapeI_shapeJ;
};

class UncollidedProblem : public LBSProblem
{
public:
  explicit UncollidedProblem(const InputParameters& params);

  ~UncollidedProblem() override;

protected:
  struct ReflectionPlane
  {
    std::uint64_t boundary_id = 0;
    Vector3 normal;
    double offset = 0.0;
  };

  struct SourcePoint
  {
    struct Subscriber
    {
      std::uint32_t cell_local_id = 0;
      double volume_weight = 1.0;
    };

    Vector3 location;
    std::vector<double> strength;
    std::vector<Subscriber> subscribers;
    std::shared_ptr<LogicalVolume> near_source_logvol;
    bool report_diagnostics = true;
    bool reflected_image = false;
    std::vector<ReflectionPlane> reflection_planes;
  };

  void PrintSimHeader() override;

  void InitializeSpatialDiscretization() override;

  void ClearBoundaries() override {}

  static Vector3 ComputeOmega(const Vector3& point0, const Vector3& point1)
  {
    double norm = (point1 - point0).Norm();
    return norm == 0. ? Vector3(0., 0., 0.) : (point1 - point0).Normalized();
  }

  /**
   * Populates cell relationships and face orientations for point source calculation.
   *
   * \param point_source The point source position vector.
   * \param cell_successors Cell successors.
   */
  void PopulateCellRelationships(const Vector3& point_source,
                                 std::vector<std::set<std::pair<size_t, double>>>& cell_successors);

  void ConstructSPLS(const std::vector<std::set<std::pair<size_t, double>>>& cell_successors,
                     const SourcePoint& source_point);

  void InitializeNearSourceRegions(const InputParameters& params);

  void InitializeReflectingBoundaries(const InputParameters& params);

  void BuildSourcePoints();

  void AddReflectedSourcePoints();

  void RaytraceNearSourceRegion(const SourcePoint& source_point);

  void ProjectReflectedImageSource(const SourcePoint& source_point);

  std::vector<double> RaytraceLine(RayTracer& ray_tracer,
                                   const Cell& cell,
                                   const Vector3& qp_xyz,
                                   const SourcePoint& source_point,
                                   double tolerance = 1e-12);

  void SweepBulkRegion(const SourcePoint& source_point);

  UncollidedMatrices ComputeUncollidedIntegrals(const Cell& cell, const Vector3& pt_loc);

  void Execute();

  void UpdateBalance(const SourcePoint& source_point);

  void AccumulateMoments(const Vector3& pt_loc);

  void WriteFluxMoments(hid_t file);

  void ComputeMoment(unsigned int ell, int m, const Vector3& pt_loc);

  void FinalizeBalance(hid_t file);

  bool IsReflectingBoundary(std::uint64_t boundary_id) const
  {
    return reflecting_boundary_ids_.count(boundary_id) != 0;
  }

  /// Near source region logical volumes.
  std::vector<std::shared_ptr<LogicalVolume>> near_source_logvols_;
  std::shared_ptr<LogicalVolume> volumetric_near_source_logvol_;
  std::vector<SourcePoint> source_points_;
  std::vector<ReflectionPlane> reflection_planes_;
  std::set<std::uint64_t> reflecting_boundary_ids_;
  /// Cell face orientations for the cells in the local cell graph.
  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;
  /// Uncollided sweep-plane local subgrid.
  std::vector<size_t> spls_;
  /// Near source uncollided sweep-plane local subgrid.
  std::vector<size_t> near_spls_;
  /// Bulk region uncollided sweep-plane local subgrid.
  std::vector<size_t> bulk_spls_;

  DenseMatrix<double> G_, M_surf_, M_;
  std::vector<Vector<double>> Phi_;

  std::vector<double> destination_phi_;

  std::string uncollided_flux_file_;
  unsigned int progress_interval_ = 5;
  unsigned int volumetric_source_quadrature_order_ = 2;
  bool project_reflected_image_sources_ = false;
  unsigned int reflected_image_projection_threads_ = 1;

  unsigned int ell_max_ = 0;
  std::vector<double> flux_moment_;
  std::vector<std::string> moment_names_;
  std::vector<std::vector<double>> accumulated_moments_;

  double production_ = 0.;
  double removal_ = 0.;
  double out_flow_ = 0.;
  bool out_flow_is_partial_ = false;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<UncollidedProblem> Create(const ParameterBlock& params);
};

} // namespace opensn
