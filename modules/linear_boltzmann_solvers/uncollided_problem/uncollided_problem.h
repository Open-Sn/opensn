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
#include <string>
#include <thread>

namespace opensn
{

class UncollidedSolver;

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
  friend class UncollidedSolver;

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
    std::vector<ReflectionPlane> reflection_planes;
  };

  enum class CellRegion
  {
    NEAR_SOURCE,
    BULK
  };

  struct Moment
  {
    unsigned int ell;
    int m;
  };

  void PrintSimHeader() override;

  void InitializeSpatialDiscretization() override;

  void ClearBoundaries() override {}

  static Vector3 ComputeOmega(const Vector3& point0, const Vector3& point1)
  {
    double norm = (point1 - point0).Norm();
    return norm == 0. ? Vector3(0., 0., 0.) : (point1 - point0).Normalized();
  }

  void BuildSweepOrdering(const SourcePoint& source_point);

  void InitializeNearSourceRegions(const InputParameters& params);

  void InitializeReflectingBoundaries(const InputParameters& params);

  void BuildSourcePoints();

  void AddReflectedSourcePoints();

  void RaytraceNearSourceRegion(const SourcePoint& source_point);

  void ProjectReflectedImageSources(unsigned int progress_interval);

  std::vector<double> RaytraceLine(RayTracer& ray_tracer,
                                   const Cell& cell,
                                   const Vector3& qp_xyz,
                                   const SourcePoint& source_point,
                                   double tolerance = 1e-12);

  void SweepBulkRegion(const SourcePoint& source_point);

  UncollidedMatrices ComputeUncollidedIntegrals(const Cell& cell, const Vector3& pt_loc);

  void Execute(const std::string& file_name, unsigned int progress_interval);

  void UpdateBalance(const SourcePoint& source_point);

  void AccumulateMoments(const Vector3& pt_loc);

  void WriteFluxMoments(hid_t file);

  void FinalizeBalance(hid_t file);

  bool IsReflectingBoundary(std::uint64_t boundary_id) const
  {
    return reflecting_boundary_ids_.count(boundary_id) != 0;
  }

  /// Near source region logical volumes.
  std::vector<std::shared_ptr<LogicalVolume>> near_source_logvols_;
  std::vector<SourcePoint> source_points_;
  std::vector<SourcePoint> reflected_source_points_;
  std::vector<ReflectionPlane> reflection_planes_;
  std::set<std::uint64_t> reflecting_boundary_ids_;
  /// Cell face orientations for the cells in the local cell graph.
  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;
  /// Near source uncollided sweep-plane local subgrid.
  std::vector<size_t> near_spls_;
  /// Bulk region uncollided sweep-plane local subgrid.
  std::vector<size_t> bulk_spls_;
  std::vector<CellRegion> cell_regions_;

  std::vector<double> destination_phi_;

  unsigned int ell_max_ = 0;
  std::vector<Moment> moments_;
  std::vector<std::string> moment_names_;
  std::vector<std::vector<double>> accumulated_moments_;

  double production_ = 0.;
  double removal_ = 0.;
  double out_flow_ = 0.;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<UncollidedProblem> Create(const ParameterBlock& params);
};

} // namespace opensn
