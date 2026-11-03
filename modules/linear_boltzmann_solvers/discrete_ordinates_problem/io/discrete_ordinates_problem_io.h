// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "framework/data_types/vector3.h"
#include "hdf5.h"
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

struct UncollidedFluxData
{
  unsigned int num_groups = 0;
  unsigned int max_moment_order = 0;
  std::uint64_t global_cell_count = 0;
  double source_rate = 0.0;
  double outflow_rate = 0.0;
  std::vector<double> local_flux_moments;
};

class DiscreteOrdinatesProblemIO
{
public:
  static UncollidedFluxData ReadUncollidedFlux(const DiscreteOrdinatesProblem& do_problem,
                                               const std::string& file_name);

  static bool ReadRestartData(DiscreteOrdinatesProblem& do_problem,
                              hid_t file_id,
                              bool allow_transient_initialization_from_steady,
                              bool repartitioned);

  static bool WriteRestartData(const DiscreteOrdinatesProblem& do_problem, hid_t file_id);

  /**
   * Write an angular flux vector to a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_src Optional angular flux source vector
   * \param role Solution role. For problem-owned data, unspecified is inferred from the problem.
   * \param orientation Adjoint orientation. For problem-owned data, unspecified is inferred.
   *
   * Each groupset records its quadrature directions and weights so that a reader can check them.
   */
  static void WriteAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src =
      std::nullopt,
    FluxSolutionRole role = FluxSolutionRole::UNSPECIFIED,
    AdjointFluxOrientation orientation = AdjointFluxOrientation::UNSPECIFIED);

  /**
   * Read an angular flux vector from a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_dest Optional angular flux destination vector
   * \param expected_role Required solution role, or unspecified to accept any role.
   * \param expected_orientation Required adjoint orientation, or unspecified to accept any.
   *
   * The file must have been written with the same rank count, partition, and quadrature; files
   * that record their quadrature are checked against the problem's. Collective: a file that
   * cannot be read or does not match the problem raises an error on every rank, and the
   * destination is left unchanged.
   */
  static void ReadAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest = std::nullopt,
    FluxSolutionRole expected_role = FluxSolutionRole::UNSPECIFIED,
    AdjointFluxOrientation expected_orientation = AdjointFluxOrientation::UNSPECIFIED);

  /**
   * Surface Angular flux
   */
  using FaceCentroidKey = std::tuple<int64_t, int64_t, int64_t>;

  /**
   * Key a surface face by its centroid. The centroid is summed over the nodes in a canonical
   * (sorted) order, so the same face yields a bit-identical key from either side of an interior
   * surface and in every file written for the same mesh. Each component is
   * `llround((centroid - centroid_origin) / centroid_spacing)`.
   *
   * \param nodes Face node locations, in any order. Must not be empty.
   * \param centroid_origin Origin subtracted from the centroid.
   * \param centroid_spacing Key spacing. Must be positive and finite.
   * \param context Description of the face source, used in error messages.
   */
  static FaceCentroidKey MakeFaceCentroidKey(std::vector<Vector3> nodes,
                                             const std::array<double, 3>& centroid_origin,
                                             double centroid_spacing,
                                             const std::string& context);

  struct SurfaceMap
  {
    std::vector<uint64_t> cell_ids;
    std::vector<uint64_t> num_face_nodes;
    /// Origin subtracted from face centroids
    std::array<double, 3> centroid_origin = {0.0, 0.0, 0.0};
    /// Spacing used to key faces by centroid: 1e-9 times the global mesh extent.
    double centroid_spacing = 1.0e-6;
    /**
     * Surface face index keyed by the face centroid, with centroid_origin subtracted from each
     * coordinate before division by centroid_spacing and rounding to the nearest integer. The
     * same face yields the same key from either side of an interior surface and in every file
     * written for the same mesh.
     */
    std::map<FaceCentroidKey, uint64_t> cell_map;
    /** Start offset into the corresponding SurfaceData::psi array for each surface cell. */
    std::vector<uint64_t> cell_stride;
    std::vector<double> nodes_x;
    std::vector<double> nodes_y;
    std::vector<double> nodes_z;
  };

  /**
   * Flattened angular-flux data for a surface and groupset.
   *
   * The angular flux `psi` is ordered by cell, face node, direction, and group. The `omega`
   * array follows cell, face node, direction, and Cartesian component, with three components per
   * node-direction pair. The `mu`, `wt_d`, and `fe_shape` arrays follow cell, face node, and
   * direction, with one value per node-direction pair. The `mass_matrix` array is ordered by cell,
   * row face node, and column face node.
   *
   * The `node_index` and `dir_index` arrays contain start offsets into `psi` for each face node and
   * node-direction pair, respectively. All arrays, including the index arrays, are empty when the
   * requested surface has no local faces.
   *
   * `psi` holds the face-node values of the problem's stored angular flux on the cell that owns the
   * face; the two sides of an interior surface are stored under separate tags. For an adjoint
   * problem the writer requires the physical orientation that an adjoint solve leaves, so the
   * stored angular flux is the adjoint flux for the listed direction and forward and adjoint data
   * with the same indices refer to the same direction.
   */
  struct SurfaceData
  {
    std::vector<double> omega;
    std::vector<double> mu;
    std::vector<double> wt_d;
    std::vector<double> mass_matrix;
    std::vector<double> fe_shape;
    std::vector<double> psi;
    std::vector<uint64_t> node_index;
    std::vector<uint64_t> dir_index;
  };

  struct SurfaceAngularFlux
  {
    int groupset_id = 0;
    /// True when the data were written from an adjoint problem.
    bool adjoint = false;
    /// Number of files in the file set, one per writing rank.
    uint64_t num_partitions = 0;
    /// Direction and group counts of the groupset, also set when the surface has no local faces.
    uint64_t num_directions = 0;
    uint64_t num_groups = 0;
    std::string surface_name;
    /// True when the tag is one side, `<name>_u` or `<name>_d`, of an interior surface; false for a
    /// boundary, whatever its name.
    bool interior = false;
    SurfaceMap mapping;
    SurfaceData data;
  };

  /**
   * Write surface angular flux vector(s) to a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param boundary_surfs Boundary surface names
   * \param interior_surfs Interior surface definitions, keyed by name, as an axis ("x", "y", or
   * "z") and a coordinate. Each surface must lie inside the mesh, must match at least one cell
   * face, must not coincide with an exterior boundary, and must not share a plane with another
   * requested surface. Cells the surface cuts through contribute no faces. Each surface is written
   * with an `_u` tag for face normals aligned with the positive specified axis and a `_d` tag for
   * the opposite orientation; neither these tags nor the surface name may equal a requested
   * boundary name. Each tag's `mesh/<tag>` group records in its `interior` attribute whether the
   * tag is one side of an interior surface or a boundary, so boundaries may have any name.
   *
   * Collective over all ranks. Surface selections may differ by rank, but each file contains only
   * its rank's tags, and DistributeSurfaceAngularFluxes requires every requested tag in every file
   * of a set. Plane-matching tolerances are relative to the global mesh extent. Each rank writes
   * `<file_base><rank>.h5`, which records whether the problem is in adjoint mode and the number of
   * files written (`num_partitions`).
   *
   * Requires `options.save_angular_flux=true`. An adjoint problem must hold its solution in the
   * physical orientation (`IsAdjointSolutionReoriented()`), which an adjoint steady-state or
   * k-eigenvalue solve leaves. Restart data are loaded in the internal sweep orientation, so an
   * adjoint solver must be executed after a restart read. Throws `std::invalid_argument` otherwise.
   */
  static void WriteSurfaceAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    const std::vector<std::string>& boundary_surfs,
    const std::map<std::string, std::pair<std::string, double>>& interior_surfs);

  /**
   * Read surface angular fluxes from one file of a file set.
   *
   * \param file_base File name base
   * \param partition File index; the file `<file_base><partition>.h5` is read.
   * \param surfaces Surface names to read. Must be non-empty. A name that is a stored tag reads
   * that tag: a boundary name, or one side of an interior surface (`_u` for face normals aligned
   * with the positive specified axis, `_d` for the opposite orientation). Any other name must be
   * an interior surface whose tags `<name>_u` and `<name>_d` both exist and are recorded as
   * interior; both are read, `_u` first. Each tag must be present in the file, possibly without
   * faces, and the names must not resolve to the same tag twice.
   *
   * Returns one entry per groupset and resolved tag, groupset-major, with `surface_name` set to
   * the tag. The file must carry every attribute
   * written by the current WriteSurfaceAngularFluxes, including `num_partitions`, which must exceed
   * `partition`. Every dataset size is checked against the face-node counts before data are read,
   * the data must be finite, and face centroid keys must be unique within each surface. The data
   * are not checked against a problem, so any rank can read any file of a set, whatever number of
   * ranks or partitioning wrote it. Not collective. Throws `std::invalid_argument` for an empty
   * `surfaces` and `std::logic_error` for a missing, malformed, or inconsistent file.
   */
  static std::vector<SurfaceAngularFlux> ReadSurfaceAngularFluxes(
    const std::string& file_base, uint64_t partition, const std::vector<std::string>& surfaces);
};

} // namespace opensn
