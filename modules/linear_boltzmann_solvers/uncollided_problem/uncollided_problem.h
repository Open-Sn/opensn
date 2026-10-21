// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/raytrace/raytracer.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector.h"
#include "framework/data_types/vector3.h"
#include "hdf5.h"
#include <atomic>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <unordered_map>

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
  ~UncollidedProblem() override;

  void ClearBoundaries() override {}

  /// Generation steps run by UncollidedSolver: locate the point sources, add their images in the
  /// reflecting planes, and compute and write the uncollided flux moments.
  void BuildSourcePoints();
  void AddReflectedSourcePoints();
  void Execute(const std::string& file_name, unsigned int progress_interval);

protected:
  /// Factory-only constructor.
  explicit UncollidedProblem(const InputParameters& params);

  /// Uncollided-generation rules.
  void CheckConfigurationErrors(std::vector<std::string>& errors) const override;

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

  static Vector3 ComputeOmega(const Vector3& point0, const Vector3& point1)
  {
    double norm = (point1 - point0).Norm();
    return norm == 0. ? Vector3(0., 0., 0.) : (point1 - point0).Normalized();
  }

  void BuildSweepOrdering(const SourcePoint& source_point);

  /// Returns the planes of the reflecting boundaries and appends an error for each invalid
  /// boundary condition.
  std::vector<ReflectionPlane> ComputeReflectionPlanes(std::vector<std::string>& errors) const;

  /// Stores the reflection planes once the configuration has been validated.
  void BuildRuntimeData() override;

  void RaytraceNearSourceRegion(const SourcePoint& source_point);

  void ProjectReflectedImageSources(unsigned int progress_interval);

  std::vector<double> RaytraceLine(RayTracer& ray_tracer,
                                   const Cell& cell,
                                   const Vector3& qp_xyz,
                                   const SourcePoint& source_point,
                                   double tolerance = 1.0e-12);

  /// Zero-allocation variant: writes into pre-allocated @p phi_out and reuses the three
  /// scratch vectors. All four buffers must be sized to @p num_groups_ before the first call
  /// (they are reset internally each call, so no manual clearing is needed between calls).
  /// On return, @p scratch_mfp holds the optical depth from the source to @p qp_xyz per group.
  void RaytraceLineInto(RayTracer& ray_tracer,
                        const Cell& cell,
                        const Vector3& qp_xyz,
                        const SourcePoint& source_point,
                        std::vector<double>& phi_out,
                        std::vector<std::pair<size_t, double>>& scratch_segs,
                        std::vector<double>& scratch_bp,
                        std::vector<double>& scratch_mfp,
                        double tolerance = 1.0e-12);

  /// Integrals of the uncollided flux of a point source over a cell.
  struct ExactCellIntegrals
  {
    /// Volume integral of the scalar flux over the cell, per group.
    std::vector<double> volume_integral;
    /// Outgoing leakage through each cell face, indexed [face][group]; zero for other faces.
    std::vector<std::vector<double>> face_leakage;
  };

  /**
   * Integrates the uncollided flux of a point source over a cell, and its leakage through each
   * outgoing face, exactly over direction.
   *
   * Every direction from the source through the cell leaves it through an outgoing face. Along
   * it the flux is S exp(-tau(r)) / (4 pi r^2), so with tau_in the optical depth from the source
   * to where the ray enters the cell and Delta the chord through the homogeneous cell, the leakage
   * integrand is S exp(-tau_in - sigma_t Delta) / (4 pi) and the volume integrand is
   * S exp(-tau_in) (1 - exp(-sigma_t Delta)) / (4 pi sigma_t). The optical depth comes from ray
   * tracing, so heterogeneous media are handled. For the cell containing the source, the ray
   * starts in the cell, so tau_in = 0 and Delta is the distance to the face.
   *
   * The directions through each PWL side triangle of each outgoing face are integrated with a
   * composite Gauss-Legendre rule in variables that keep the integrands smooth however close the
   * source is to a face, edge, or vertex, at points that lie on the face.
   *
   * Returns std::nullopt when a ray from the source could cross the cell more than once: for a
   * cell that does not contain the source, when it is not convex; and for the cell containing
   * the source, when a face is not directly visible from it. Also returns std::nullopt when the
   * source is on the inner side of a face that is not outgoing.
   */
  std::optional<ExactCellIntegrals>
  IntegrateCellExactly(RayTracer& ray_tracer, const Cell& cell, const SourcePoint& source_point);

  /**
   * Returns the fraction of the full solid angle that a cell subtends at a point on or in it.
   *
   * The signed solid angles of the PWL side triangles of the cell's faces are summed, so the
   * result is exact for cells of any shape, including non-convex ones. It is 1 for a point
   * inside the cell, 1/2 on a planar face, and smaller on an edge or vertex.
   */
  double SubtendedAngleFraction(const Cell& cell, const Vector3& point) const;

  /// True when every vertex of the cell lies on the inner side of the plane of every face.
  bool IsConvex(const Cell& cell) const;

  void SweepBulkRegion(const SourcePoint& source_point);

  UncollidedMatrices ComputeUncollidedIntegrals(const Cell& cell, const Vector3& pt_loc);

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
  /// Boundary conditions as (name, type).
  std::vector<std::pair<std::string, std::string>> boundary_conditions_;
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
  /// Bounding-box diagonal for each local cell, pre-computed once for RayTracer tolerance scaling.
  std::vector<double> cell_sizes_;
  /// False when the mesh contains any geometrically non-convex cell; true for standard hex/tet
  /// meshes.
  bool all_cells_convex_ = true;

  /// Precomputed face vertex data -- bypasses the std::map-backed VertexHandler in the hot path.
  /// Supports faces with up to max_sides vertices (covers common polyhedral mesh faces).
  struct FaceVertData
  {
    static constexpr uint32_t max_sides = 8;
    Vector3 verts[max_sides]; // inline vertex positions (no tree traversal)
    Vector3 centroid;
    uint64_t neighbor_id = 0;
    uint32_t num_sides = 0;
    uint32_t pad = 0;
  };
  std::vector<FaceVertData>
    all_face_verts_; // flat: [cell_0_face_0, cell_0_face_1, ..., cell_1_face_0, ...]
  std::vector<uint32_t> cell_face_offsets_; // [num_cells] first face index in all_face_verts_
  std::vector<uint32_t> cell_num_faces_;    // [num_cells] number of faces per cell
  /// O(1) global-to-local cell-ID lookup; replaces two O(log N) std::map calls per while-loop step.
  std::unordered_map<uint64_t, uint32_t> global_to_local_id_;
  /// True when all_cells_convex_ and all faces have <= FaceVertData::max_sides vertices.
  bool use_fast_trace_ = false;

  unsigned int ell_max_ = 0;
  /// Integrate the leakages and volume integrals of the cells that share a vertex with a point
  /// source's cell exactly over direction instead of with spatial quadrature.
  bool exact_near_source_leakage_ = false;
  std::vector<Moment> moments_;
  std::vector<std::string> moment_names_;
  std::vector<std::vector<double>> accumulated_moments_;

  double production_ = 0.;
  double physical_removal_ = 0.;
  double reflected_removal_ = 0.;
  double out_flow_ = 0.;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<UncollidedProblem> Create(const ParameterBlock& params);
};

} // namespace opensn
