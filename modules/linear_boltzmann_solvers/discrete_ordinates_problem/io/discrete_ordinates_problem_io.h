// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

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
                              bool allow_transient_initialization_from_steady);

  static bool WriteRestartData(const DiscreteOrdinatesProblem& do_problem, hid_t file_id);

  /**
   * Write an angular flux vector to a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_src Optional angular flux source vector
   */
  static void WriteAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src =
      std::nullopt);

  /**
   * Read an angular flux vector from a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_dest Optional angular flux destination vector
   */
  static void ReadAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest =
      std::nullopt);

  /**
   * Surface Angular flux
   */
  using FaceCentroidKey = std::tuple<int64_t, int64_t, int64_t>;

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
   * face; the two sides of an interior surface are stored under separate tags. After an adjoint
   * steady-state solve the stored angular flux is the adjoint flux for the listed direction, so
   * forward and adjoint data with the same indices refer to the same direction.
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
    std::string surface_name;
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
   * "z") and a coordinate. Each surface must match at least one cell face, must lie on cell faces
   * inside the mesh, must not coincide with an exterior boundary, and must not share a plane with
   * another requested surface. Each surface is written with an `_u` tag for face normals aligned
   * with the positive specified axis and a `_d` tag for the opposite orientation; these tags must
   * not equal a requested boundary name.
   *
   * Collective over all ranks. Surface selections may differ by rank. Plane-matching tolerances
   * are relative to the global mesh extent. Each rank writes `<file_base><rank>.h5`, which records
   * whether the problem is in adjoint mode.
   */
  static void WriteSurfaceAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    const std::vector<std::string>& boundary_surfs,
    const std::map<std::string, std::pair<std::string, double>>& interior_surfs);

  /**
   * Read a surface angular flux vector from a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param surfaces Stored surface tags to read. Boundary tags are their boundary names. Interior
   * surface tags use `_u` for face normals aligned with the positive specified axis and `_d` for
   * the opposite orientation.
   *
   * Reads `<file_base><rank>.h5` on each rank; the files must have been written with the same mesh
   * partitioning. Not collective.
   */
  static std::vector<SurfaceAngularFlux>
  ReadSurfaceAngularFluxes(DiscreteOrdinatesProblem& do_problem,
                           const std::string& file_base,
                           const std::vector<std::string>& surfaces);
};

} // namespace opensn
