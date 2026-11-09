// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "framework/runtime.h"
#include <array>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace opensn
{

/// Default per-rank send budget, in bytes, of one surface-data exchange.
const std::uint64_t surface_redistribution_bytes = std::uint64_t{64} << 20;

/// Quadrature of one groupset of a surface angular-flux file set.
struct SurfaceQuadrature
{
  std::uint64_t num_dirs = 0;
  std::uint64_t num_groups = 0;
  /// First group of the groupset.
  std::uint64_t first_group = 0;
  /// Direction vectors, three Cartesian components per direction.
  std::vector<double> omegas;
  /// Quadrature weights, one per direction.
  std::vector<double> weights;
};

/**
 * The faces of one stored surface tag that are owned by this rank.
 *
 * A face is owned by rank `SurfaceFaceOwner(key, P)`, where `key` is its centroid key, so every
 * data set read on the same communicator places a given physical face, from either side of an
 * interior surface, on the same rank. Faces are stored in arrival order; `index` gives the local
 * face index of each key in key order.
 */
struct SurfaceShard
{
  struct GroupsetValues
  {
    /// Direction cosine with the face normal, ordered by face and direction.
    std::vector<double> mu;
    /// Angular flux, ordered by face, face node, direction, and group. For an interior tag, the
    /// values of the face's cell; for a boundary tag, the upwind values written by
    /// WriteSurfaceAngularFluxes (`upwind_flux`), which take the boundary condition for the
    /// directions that enter the domain.
    std::vector<double> psi;
  };

  /// True when the tag is one side, `<name>_u` or `<name>_d`, of an interior surface; false for a
  /// boundary, whatever its name. Recorded by the writer and identical on all ranks.
  bool interior = false;
  /// (face key, local face index) pairs sorted by key.
  std::vector<std::pair<DiscreteOrdinatesProblemIO::FaceCentroidKey, std::uint64_t>> index;
  /// Start offset of each face's nodes. Its size is the number of faces plus one.
  std::vector<std::uint64_t> node_offset = {0};
  /// Start offset of each face's mass matrix. Its size is the number of faces plus one.
  std::vector<std::uint64_t> mass_offset = {0};
  /// Face node coordinates, three per node.
  std::vector<double> coords;
  /// Face mass matrices, row-major by face node, including the curvilinear measure scale.
  std::vector<double> mass;
  /// Values per groupset.
  std::vector<GroupsetValues> groupsets;

  std::uint64_t GetNumFaces() const { return node_offset.size() - 1; }

  std::uint64_t GetNumFaceNodes(const std::uint64_t face) const
  {
    return node_offset[face + 1] - node_offset[face];
  }

  /// Return the local index of the face with `key`, if this rank owns it.
  std::optional<std::uint64_t> Find(const DiscreteOrdinatesProblemIO::FaceCentroidKey& key) const;
};

/// Surface angular-flux data distributed by face across the ranks of a communicator.
struct DistributedSurfaceData
{
  /// True when the data were written from an adjoint problem.
  bool adjoint = false;
  /// Number of files the data were read from.
  std::uint64_t num_partitions = 0;
  std::array<double, 3> centroid_origin = {0.0, 0.0, 0.0};
  double centroid_spacing = 0.0;
  /// Quadrature per groupset. Identical on all ranks.
  std::vector<SurfaceQuadrature> groupsets;
  /// Locally owned faces keyed by stored surface tag. Every rank holds every tag.
  std::map<std::string, SurfaceShard> surfaces;
};

/// How a surface of the forward data is evaluated; see ComputeSurfaceResponses.
enum class SurfaceKind
{
  /// An exterior boundary, stored under its boundary name.
  boundary,
  /// An interior surface, stored as its two sides `<name>_u` and `<name>_d`.
  interior
};

/// The surface response of one boundary or interior surface.
struct SurfaceResponse
{
  /// The surface as requested: a boundary name, or an interior surface name, bare (`X`) or with
  /// the side tag that sets its normal (`X_u` or `X_d`).
  std::string name;
  SurfaceKind kind = SurfaceKind::boundary;
  /// Contribution of the directions with positive cosine with the surface normal. For an interior
  /// surface, the normal is the outward normal of the selected side, `X_u` (the positive axis) for
  /// `X` and `X_u` and `X_d` (the negative axis) for `X_d`, and these directions leave the cell of
  /// that side. For a boundary, the normal is the outward normal and these are the outflow
  /// directions.
  double outgoing = 0.0;
  /// Contribution of the directions with negative cosine with the normal. For an interior surface,
  /// these leave the cell of the other side and cross back. For a boundary, they are the inflow
  /// directions: the forward boundary-condition value (a forward boundary source; 0 for a vacuum
  /// boundary) against the adjoint cell value.
  double incoming = 0.0;

  /// The response carried across the surface along its normal: `outgoing + incoming`.
  double Total() const { return outgoing + incoming; }
};

/// Rank that owns the face with `key` on a communicator with `num_ranks` ranks.
int SurfaceFaceOwner(const DiscreteOrdinatesProblemIO::FaceCentroidKey& key, int num_ranks);

/**
 * Distribute the faces of a surface angular-flux file set written by
 * `DiscreteOrdinatesProblemIO::WriteSurfaceAngularFluxes` across the ranks of `comm`.
 *
 * Each file is read with `DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes`, which does not
 * check the data against a problem. The files `<file_base>0.h5`, ...,
 * `<file_base>(N-1).h5` are found from the `num_partitions` attribute of the first file. N need not
 * equal the size of `comm`, and the data need not match the partitioning of any current mesh. File
 * k is read by rank k modulo the size of `comm`, one file per rank at a time, and its faces are
 * sent to their owner ranks (see SurfaceFaceOwner) in exchanges in which each rank sends at most
 * `max_send_bytes` packed bytes (at least one face), before the next file is read. The peak memory
 * of a rank is therefore one file's requested surfaces, as returned by the reader, plus its owned
 * faces and one exchange. A rank receives from every rank in an exchange, so its receive buffer is
 * bounded by the number of ranks times `max_send_bytes`; faces spread evenly over owners, so it is
 * typically about `max_send_bytes`. Every file must contain every requested tag, possibly
 * without faces.
 *
 * \param file_base File name base.
 * \param surfaces Surface names to read. Must be non-empty and unique. A boundary name reads that
 *        tag, and the name of an interior surface reads both of its tags, `<name>_u` and
 *        `<name>_d`. A single side's tag is rejected unless its other side is also requested,
 *        because surface responses need both sides. Names are resolved against file 0, and every
 *        other file is read with the resolved tags.
 * \param comm Communicator to distribute the faces over.
 * \param max_send_bytes Per-rank send budget of one exchange.
 *
 * Collective over `comm`. Throws `std::invalid_argument` for invalid arguments and
 * `std::logic_error`, on all ranks, for missing, inconsistent, or malformed files, for a tag with
 * no faces in any file, and for duplicate faces.
 */
DistributedSurfaceData
DistributeSurfaceAngularFluxes(const std::string& file_base,
                               const std::vector<std::string>& surfaces,
                               const mpi::Communicator& comm,
                               std::uint64_t max_send_bytes = surface_redistribution_bytes);

/**
 * Evaluate the surface inner products of forward and adjoint angular fluxes.
 *
 * For a surface S with faces f,
 *   R = sum_{f in S} sum_{d in D_f} mu_{f,d} w_d sum_ij M_{f,ij} sum_g
 *         psi_f(i, d, g) psi^dagger_f(j, d, g),
 * with mu_{f,d} the cosine of direction d with the normal n_f of face f. The forward tags of the
 * forward data define the surfaces, classified by the kind the writer records for each tag
 * (SurfaceShard::interior), not by its name:
 * - A boundary: n_f is the outward normal, D_f holds every direction, and both fluxes are upwind
 *   face values from the boundary tag of each data set. For the outflow directions (mu_{f,d} > 0)
 *   they are the forward flux of the face's cell and the adjoint boundary-condition value (0 for a
 *   vacuum boundary), which is upwind for the adjoint because it travels along -Omega_d; for the
 *   inflow directions, the forward boundary-condition value and the adjoint flux of the cell. This
 *   is the boundary term of the discrete adjoint identity, so it is exact for the discretization
 *   like an interior surface, and the totals of all boundaries sum to <q, psi^dagger> -
 *   <q^dagger, psi> for the volumetric sources q and q^dagger.
 * - An interior surface, with both interior tags `X_u` and `X_d` in the forward data: the faces are
 * those of `X_u`, n_f is the `X_u` normal (the positive axis), and D_f holds every direction. The
 * forward flux comes from the side upwind of the direction (`X_u` for mu_{f,d} > 0, `X_d` for
 *   mu_{f,d} < 0) and the adjoint flux from the other side, which is upwind for the adjoint.
 *   Requested as `X` or `X_u`, R is measured along the `X_u` normal; requested as `X_d`, along the
 *   `X_d` normal (the negative axis), which gives -R and swaps and negates `outgoing` and
 *   `incoming`.
 * Each response is split into `outgoing` (mu_{f,d} > 0) and `incoming` (mu_{f,d} < 0) parts; the
 * response is their sum (SurfaceResponse::Total). OpenSn quadrature weights sum to 1, so for an
 * adjoint problem driven by a volumetric, material, or point response source, R equals the response
 * it carries across the surface. Forward and adjoint directions are paired by equal direction
 * vectors and must have equal weights. Face nodes are paired by location.
 *
 * \param forward Forward data. Must not be adjoint.
 * \param adjoint Adjoint data. Must be adjoint and must be distributed over `comm`.
 * \param surfaces Surfaces to evaluate: boundary names and interior surface names (`X`, or `X_u`
 *        or `X_d` to choose the normal). Each surface may be requested once. If empty, every
 * surface of the forward data is evaluated, interior surfaces along their `X_u` normals.
 * \param comm Communicator both data sets were distributed over.
 * \return One surface response per surface, in name order when `surfaces` is empty and in the given
 * order otherwise. Values are global sums.
 *
 * Collective over `comm`. Throws `std::invalid_argument` for incompatible data sets, unknown
 * surfaces, a surface requested more than once, a side of an interior surface without its other
 * side, or an adjoint tag whose kind differs from its
 * forward surface, and `std::logic_error`, on all ranks, if any
 * face lacks a matching face in another required tag or matching faces have different node
 * locations.
 */
std::vector<SurfaceResponse> ComputeSurfaceResponses(const DistributedSurfaceData& forward,
                                                     const DistributedSurfaceData& adjoint,
                                                     const std::vector<std::string>& surfaces,
                                                     const mpi::Communicator& comm);

} // namespace opensn
