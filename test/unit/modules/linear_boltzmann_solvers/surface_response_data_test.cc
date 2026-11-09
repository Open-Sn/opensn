// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/response_evaluator/surface_response_data.h"
#include "framework/runtime.h"
#include "framework/utils/hdf_utils.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

using namespace opensn;

namespace
{

/*
 * Synthetic surface angular-flux files in the WriteSurfaceAngularFluxes format.
 *
 * Each surface tag lies on a plane x = const and consists of num_faces two-node segments of length
 * face_length along y, so the surface spans y in [0, surface_length]. The angular flux is linear in
 * y, so the face mass matrix face_length/6 [[2, 1], [1, 2]] integrates the product of two fluxes
 * exactly and the response has the closed form of ExpectedResponse. Boundary tags also store upwind
 * values, which differ from the cell values ('i' instead of the tag's side) for the directions
 * that enter the domain, so the boundary response checks that the adjoint upwind values are used.
 */
const int num_faces = 7;
const double face_length = 0.2;
const double surface_length = num_faces * face_length;
const std::vector<std::array<double, 3>> directions = {
  {0.6, 0.8, 0.0}, {-0.6, 0.8, 0.0}, {-0.6, -0.8, 0.0}, {0.6, -0.8, 0.0}};
const std::vector<double> weights = {0.1, 0.2, 0.3, 0.4};
/// Groups per groupset.
const std::vector<std::uint64_t> groups_per_groupset = {2, 1};

/// Angular flux a + b y of a data set (`adjoint`) on one side ('u', 'd', or 'b') of a surface, or
/// the boundary-condition value ('i') of a boundary.
std::array<double, 2>
FluxCoefficients(const bool adjoint,
                 const char side,
                 const std::array<double, 3>& omega,
                 const size_t gs,
                 const std::uint64_t g)
{
  const double s = side == 'u' ? 1.0 : (side == 'd' ? 2.0 : (side == 'b' ? 3.0 : 4.0));
  const double a = (adjoint ? 0.5 : 1.0) + 0.3 * s + omega[0] + 0.7 * omega[1] +
                   0.1 * static_cast<double>(g) + 0.05 * static_cast<double>(gs);
  const double b = (adjoint ? -0.4 : 0.6) + 0.2 * s * omega[1] + 0.1 * static_cast<double>(g);
  return {a, b};
}

/// Integral over [0, surface_length] of (a1 + b1 y)(a2 + b2 y).
double
IntegrateProduct(const std::array<double, 2>& f, const std::array<double, 2>& g)
{
  return f[0] * g[0] * surface_length +
         (f[0] * g[1] + f[1] * g[0]) * surface_length * surface_length / 2.0 +
         f[1] * g[1] * surface_length * surface_length * surface_length / 3.0;
}

/**
 * Closed-form response of an interior surface with normal +x (`interior`) or a boundary with
 * outward normal +x. On the boundary each flux takes its upwind value: the cell value for the
 * directions leaving the domain along its transport direction and the boundary-condition value for
 * those entering it.
 */
double
ExpectedResponse(const bool interior)
{
  double response = 0.0;
  for (size_t gs = 0; gs < groups_per_groupset.size(); ++gs)
    for (size_t d = 0; d < directions.size(); ++d)
    {
      const double mu = directions[d][0];
      for (std::uint64_t g = 0; g < groups_per_groupset[gs]; ++g)
      {
        if (not interior)
          // Upwind values: the forward cell value and adjoint boundary value for mu > 0, the
          // forward boundary value and adjoint cell value for mu < 0.
          response +=
            mu * weights[d] *
            IntegrateProduct(FluxCoefficients(false, mu > 0.0 ? 'b' : 'i', directions[d], gs, g),
                             FluxCoefficients(true, mu > 0.0 ? 'i' : 'b', directions[d], gs, g));
        else if (mu > 0.0)
          response += mu * weights[d] *
                      IntegrateProduct(FluxCoefficients(false, 'u', directions[d], gs, g),
                                       FluxCoefficients(true, 'd', directions[d], gs, g));
        else
          response += mu * weights[d] *
                      IntegrateProduct(FluxCoefficients(false, 'd', directions[d], gs, g),
                                       FluxCoefficients(true, 'u', directions[d], gs, g));
      }
    }
  return response;
}

/// Closed-form incoming (mu < 0) part of the boundary response: the forward boundary value against
/// the adjoint cell value.
double
ExpectedBoundaryIncoming()
{
  double incoming = 0.0;
  for (size_t gs = 0; gs < groups_per_groupset.size(); ++gs)
    for (size_t d = 0; d < directions.size(); ++d)
      if (directions[d][0] < 0.0)
        for (std::uint64_t g = 0; g < groups_per_groupset[gs]; ++g)
          incoming += directions[d][0] * weights[d] *
                      IntegrateProduct(FluxCoefficients(false, 'i', directions[d], gs, g),
                                       FluxCoefficients(true, 'b', directions[d], gs, g));
  return incoming;
}

/// Closed-form outgoing (mu > 0) and incoming (mu < 0) parts of the interior surface response
/// along +x.
std::pair<double, double>
ExpectedInteriorParts()
{
  double outgoing = 0.0;
  double incoming = 0.0;
  for (size_t gs = 0; gs < groups_per_groupset.size(); ++gs)
    for (size_t d = 0; d < directions.size(); ++d)
    {
      const double mu = directions[d][0];
      for (std::uint64_t g = 0; g < groups_per_groupset[gs]; ++g)
        (mu > 0.0 ? outgoing : incoming) +=
          mu * weights[d] *
          IntegrateProduct(FluxCoefficients(false, mu > 0.0 ? 'u' : 'd', directions[d], gs, g),
                           FluxCoefficients(true, mu > 0.0 ? 'd' : 'u', directions[d], gs, g));
    }
  return {outgoing, incoming};
}

struct TagSpec
{
  std::string name;
  double x = 0.0;
  /// x-component of the stored face normal.
  double normal_x = 1.0;
  char side = 'u';
  /// Faces written to each file.
  std::vector<std::vector<int>> faces_per_file;
  /// Surface kind recorded in the files; the default is an interior side unless `side` is 'b'.
  std::optional<bool> recorded_interior;
};

struct SetSpec
{
  bool adjoint = false;
  /// Omit attributes that earlier versions of OpenSn did not write.
  bool write_num_partitions = true;
  bool write_adjoint = true;
  bool write_interior = true;
  bool write_file_set_id = true;
  bool write_first_group = true;
  bool write_upwind = true;
  /// Identifier of the set, and an optional (file, identifier) pair for a file of another set.
  std::uint64_t file_set_id = 7;
  std::optional<std::pair<size_t, std::uint64_t>> stray_file;
  /// Added to the first group of every groupset.
  std::uint64_t first_group_shift = 0;
  double centroid_spacing = 1.0e-9;
  /// Store each face's nodes in reverse order.
  bool reverse_nodes = false;
  /// Store directions in reverse order.
  bool reverse_directions = false;
  /// Factor applied to the fixture quadrature weights.
  double weight_scale = 1.0;
  /// A file whose quadrature weights are doubled, unlike the rest of the set.
  std::optional<size_t> odd_quadrature_file;
  /// Drop the last angular-flux value of every non-empty tag, so its size is inconsistent.
  bool truncate_flux = false;
  /// Store a NaN as the first angular-flux value of every non-empty tag.
  bool nan_flux = false;
  /// Number of files recorded in each file; the default is the number written.
  std::optional<std::uint64_t> recorded_partitions;
  std::vector<TagSpec> tags;
};

void
WriteSet(const std::string& file_base, const SetSpec& spec)
{
  const auto num_files = spec.tags.front().faces_per_file.size();
  for (size_t file_id = 0; file_id < num_files; ++file_id)
  {
    const auto file_name = file_base + std::to_string(file_id) + ".h5";
    const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(file.Id(), 0);
    const auto id = file.Id();
    ASSERT_TRUE(H5CreateAttribute(id, "num_groupsets", std::uint64_t{groups_per_groupset.size()}));
    if (spec.write_adjoint)
    {
      ASSERT_TRUE(H5CreateAttribute<bool>(id, "adjoint", spec.adjoint));
    }
    if (spec.write_file_set_id)
    {
      const auto set_id = spec.stray_file.has_value() and spec.stray_file->first == file_id
                            ? spec.stray_file->second
                            : spec.file_set_id;
      ASSERT_TRUE(H5CreateAttribute(id, "file_set_id", set_id));
    }
    if (spec.write_num_partitions)
    {
      ASSERT_TRUE(H5CreateAttribute(
        id, "num_partitions", spec.recorded_partitions.value_or(std::uint64_t{num_files})));
    }
    for (size_t d = 0; d < 3; ++d)
      ASSERT_TRUE(H5CreateAttribute(id, "centroid_origin_" + std::to_string(d), 0.0));
    ASSERT_TRUE(H5CreateAttribute(id, "centroid_spacing", spec.centroid_spacing));
    ASSERT_TRUE(H5CreateGroup(id, "mesh"));

    for (const auto& tag : spec.tags)
    {
      const auto& faces = tag.faces_per_file[file_id];
      std::vector<std::uint64_t> cell_ids, num_face_nodes;
      std::vector<double> x, y, z, normals, fe_shape, mass;
      for (const auto face : faces)
      {
        cell_ids.push_back(static_cast<std::uint64_t>(face));
        num_face_nodes.push_back(2);
        normals.insert(normals.end(), {tag.normal_x, 0.0, 0.0});
        fe_shape.insert(fe_shape.end(), {face_length / 2.0, face_length / 2.0});
        mass.insert(mass.end(),
                    {face_length / 3.0, face_length / 6.0, face_length / 6.0, face_length / 3.0});
        std::array<double, 2> ys = {face * face_length, (face + 1) * face_length};
        if (spec.reverse_nodes)
          std::swap(ys[0], ys[1]);
        for (const auto value : ys)
        {
          x.push_back(tag.x);
          y.push_back(value);
          z.push_back(0.0);
        }
      }
      const auto mesh = "mesh/" + tag.name;
      ASSERT_TRUE(H5CreateGroup(id, mesh));
      if (spec.write_interior)
      {
        ASSERT_TRUE(H5CreateAttribute<bool>(
          id, mesh + "/interior", tag.recorded_interior.value_or(tag.side != 'b')));
      }
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/cell_ids", cell_ids));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/num_face_nodes", num_face_nodes));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/nodes_x", x));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/nodes_y", y));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/nodes_z", z));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/normals", normals));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/fe_shape", fe_shape));
      ASSERT_TRUE(H5WriteDataset1D(id, mesh + "/M_ij", mass));
    }

    for (size_t gs = 0; gs < groups_per_groupset.size(); ++gs)
    {
      const auto group = "groupset_" + std::to_string(gs);
      ASSERT_TRUE(H5CreateGroup(id, group));
      ASSERT_TRUE(
        H5CreateAttribute(id, group + "/num_directions", std::uint64_t{directions.size()}));
      ASSERT_TRUE(H5CreateAttribute(id, group + "/num_groups", groups_per_groupset[gs]));
      if (spec.write_first_group)
      {
        std::uint64_t first_group = spec.first_group_shift;
        for (size_t prior = 0; prior < gs; ++prior)
          first_group += groups_per_groupset[prior];
        ASSERT_TRUE(H5CreateAttribute(id, group + "/first_group", first_group));
      }
      std::vector<double> omega, wt;
      for (size_t k = 0; k < directions.size(); ++k)
      {
        const auto d = spec.reverse_directions ? directions.size() - 1 - k : k;
        omega.insert(omega.end(), directions[d].begin(), directions[d].end());
        const double odd = spec.odd_quadrature_file == file_id ? 2.0 : 1.0;
        wt.push_back(odd * spec.weight_scale * weights[d]);
      }
      ASSERT_TRUE(H5WriteDataset1D(id, group + "/omega", omega));
      ASSERT_TRUE(H5WriteDataset1D(id, group + "/weights", wt));
      for (const auto& tag : spec.tags)
      {
        const bool boundary = not tag.recorded_interior.value_or(tag.side != 'b');
        std::vector<double> psi, upwind;
        for (const auto face : tag.faces_per_file[file_id])
        {
          std::array<double, 2> ys = {face * face_length, (face + 1) * face_length};
          if (spec.reverse_nodes)
            std::swap(ys[0], ys[1]);
          for (const auto node_y : ys)
            for (size_t k = 0; k < directions.size(); ++k)
            {
              const auto d = spec.reverse_directions ? directions.size() - 1 - k : k;
              const double mu_d = directions[d][0] * tag.normal_x;
              const bool enters = spec.adjoint ? mu_d > 0.0 : mu_d < 0.0;
              for (std::uint64_t g = 0; g < groups_per_groupset[gs]; ++g)
              {
                const auto c = FluxCoefficients(spec.adjoint, tag.side, directions[d], gs, g);
                psi.push_back(c[0] + c[1] * node_y);
                const auto u =
                  enters ? FluxCoefficients(spec.adjoint, 'i', directions[d], gs, g) : c;
                upwind.push_back(u[0] + u[1] * node_y);
              }
            }
        }
        if (spec.truncate_flux and not psi.empty())
          psi.pop_back();
        if (spec.nan_flux and not psi.empty())
          psi.front() = std::numeric_limits<double>::quiet_NaN();
        const auto surf = group + "/" + tag.name;
        ASSERT_TRUE(H5CreateGroup(id, surf));
        ASSERT_TRUE(H5WriteDataset1D(id, surf + "/surf_flux", psi));
        if (boundary and spec.write_upwind)
        {
          if (spec.truncate_flux and not upwind.empty())
            upwind.pop_back();
          ASSERT_TRUE(H5WriteDataset1D(id, surf + "/upwind_flux", upwind));
        }
      }
    }
  }
}

/// Forward set: five files, faces split unevenly, the two sides of the plane in different files.
SetSpec
ForwardSpec()
{
  SetSpec spec;
  spec.tags = {
    {"plane_u", 0.5, 1.0, 'u', {{0, 1, 2, 3, 4}, {}, {5}, {}, {6}}},
    {"plane_d", 0.5, -1.0, 'd', {{6}, {0, 1}, {}, {2, 3, 4, 5}, {}}},
    {"right", 1.0, 1.0, 'b', {{}, {0, 1, 2, 3, 4, 5, 6}, {}, {}, {}}},
  };
  return spec;
}

/// Adjoint set: three files with another split, reversed node order, and reordered directions.
SetSpec
AdjointSpec()
{
  SetSpec spec;
  spec.adjoint = true;
  spec.reverse_nodes = true;
  spec.reverse_directions = true;
  spec.tags = {
    {"plane_u", 0.5, 1.0, 'u', {{3, 4, 5, 6}, {0}, {1, 2}}},
    {"plane_d", 0.5, -1.0, 'd', {{}, {0, 1, 2, 3, 4, 5, 6}, {}}},
    {"right", 1.0, 1.0, 'b', {{5, 6}, {0, 1, 2}, {3, 4}}},
  };
  return spec;
}

class SurfaceResponseDataTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    // All ranks use rank 0's process id so concurrent test runs do not share files.
    auto pid = static_cast<std::uint64_t>(getpid());
    mpi_comm.broadcast(pid, 0);
    directory_ =
      std::filesystem::temp_directory_path() / ("opensn_surface_response_" + std::to_string(pid));
    if (mpi_comm.rank() == 0)
      std::filesystem::create_directories(directory_);
    mpi_comm.barrier();
  }

  void TearDown() override
  {
    mpi_comm.barrier();
    if (mpi_comm.rank() == 0)
      std::filesystem::remove_all(directory_);
  }

  /// Write a set on rank 0 and return its file base.
  std::string Write(const std::string& name, const SetSpec& spec)
  {
    const auto file_base = (directory_ / (name + "_p")).string();
    if (mpi_comm.rank() == 0)
      WriteSet(file_base, spec);
    mpi_comm.barrier();
    return file_base;
  }

  static std::uint64_t GlobalFaces(const DistributedSurfaceData& data, const std::string& tag)
  {
    auto faces = data.surfaces.at(tag).GetNumFaces();
    mpi_comm.all_reduce(faces, mpi::op::sum<std::uint64_t>());
    return faces;
  }

  std::filesystem::path directory_;
};

const std::vector<std::string> all_tags = {"plane_u", "plane_d", "right"};

} // namespace

TEST_F(SurfaceResponseDataTest, MatchesClosedFormAcrossPartitions)
{
  const auto forward_base = Write("forward", ForwardSpec());
  const auto adjoint_base = Write("adjoint", AdjointSpec());

  const auto forward = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(adjoint_base, all_tags, mpi_comm);
  EXPECT_EQ(forward.num_partitions, 5U);
  EXPECT_EQ(adjoint.num_partitions, 3U);
  EXPECT_FALSE(forward.adjoint);
  EXPECT_TRUE(adjoint.adjoint);
  for (const auto& tag : all_tags)
  {
    EXPECT_EQ(GlobalFaces(forward, tag), static_cast<std::uint64_t>(num_faces));
    EXPECT_EQ(GlobalFaces(adjoint, tag), static_cast<std::uint64_t>(num_faces));
  }

  const auto surface_responses = ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm);
  ASSERT_EQ(surface_responses.size(), 2U);
  // Name order: "plane" (interior) before "right" (boundary).
  EXPECT_EQ(surface_responses[0].name, "plane");
  EXPECT_EQ(surface_responses[0].kind, SurfaceKind::interior);
  EXPECT_EQ(surface_responses[1].name, "right");
  EXPECT_EQ(surface_responses[1].kind, SurfaceKind::boundary);
  // The boundary's incoming part pairs the forward boundary value with the adjoint cell value.
  const double boundary_incoming = ExpectedBoundaryIncoming();
  EXPECT_LT(boundary_incoming, 0.0);
  EXPECT_NEAR(
    surface_responses[1].incoming, boundary_incoming, 1.0e-13 * std::abs(boundary_incoming));

  // Sums of a few hundred O(1) products: allow a small multiple of machine precision.
  const double interior = ExpectedResponse(true);
  const double boundary = ExpectedResponse(false);
  EXPECT_NEAR(surface_responses[0].outgoing + surface_responses[0].incoming,
              interior,
              1.0e-13 * std::abs(interior));
  EXPECT_NEAR(surface_responses[1].Total(), boundary, 1.0e-13 * std::abs(boundary));
  // Both interior directions contribute, with opposite signs of mu.
  EXPECT_GT(surface_responses[0].outgoing, 0.0);
  EXPECT_LT(surface_responses[0].incoming, 0.0);

  // Selecting one surface gives the same surface response.
  const auto selected = ComputeSurfaceResponses(forward, adjoint, {"right"}, mpi_comm);
  ASSERT_EQ(selected.size(), 1U);
  EXPECT_DOUBLE_EQ(selected[0].outgoing, surface_responses[1].outgoing);
}

TEST_F(SurfaceResponseDataTest, SmallSendBudgetGivesSameResult)
{
  const auto forward_base = Write("forward", ForwardSpec());
  const auto adjoint_base = Write("adjoint", AdjointSpec());

  // A one-byte budget sends one face per rank per round.
  const auto forward = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm, 1);
  const auto adjoint = DistributeSurfaceAngularFluxes(adjoint_base, all_tags, mpi_comm, 1);
  const auto small = ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm);

  const auto forward_default = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm);
  const auto adjoint_default = DistributeSurfaceAngularFluxes(adjoint_base, all_tags, mpi_comm);
  const auto reference = ComputeSurfaceResponses(forward_default, adjoint_default, {}, mpi_comm);

  // Owned faces are evaluated in key order, so the chunking does not change the sum.
  ASSERT_EQ(small.size(), reference.size());
  for (size_t s = 0; s < small.size(); ++s)
  {
    EXPECT_EQ(small[s].outgoing, reference[s].outgoing);
    EXPECT_EQ(small[s].incoming, reference[s].incoming);
  }
}

TEST_F(SurfaceResponseDataTest, RejectsFilesFromEarlierVersions)
{
  auto without_partitions = ForwardSpec();
  without_partitions.write_num_partitions = false;
  EXPECT_THROW(
    DistributeSurfaceAngularFluxes(Write("no_partitions", without_partitions), all_tags, mpi_comm),
    std::logic_error);

  auto without_adjoint = AdjointSpec();
  without_adjoint.write_adjoint = false;
  EXPECT_THROW(
    DistributeSurfaceAngularFluxes(Write("no_adjoint", without_adjoint), all_tags, mpi_comm),
    std::logic_error);

  auto without_kinds = ForwardSpec();
  without_kinds.write_interior = false;
  EXPECT_THROW(DistributeSurfaceAngularFluxes(Write("no_kinds", without_kinds), all_tags, mpi_comm),
               std::logic_error);

  auto without_set_id = ForwardSpec();
  without_set_id.write_file_set_id = false;
  EXPECT_THROW(
    DistributeSurfaceAngularFluxes(Write("no_set_id", without_set_id), all_tags, mpi_comm),
    std::logic_error);

  auto without_first_group = ForwardSpec();
  without_first_group.write_first_group = false;
  EXPECT_THROW(DistributeSurfaceAngularFluxes(
                 Write("no_first_group", without_first_group), all_tags, mpi_comm),
               std::logic_error);

  auto without_upwind = AdjointSpec();
  without_upwind.write_upwind = false;
  EXPECT_THROW(
    DistributeSurfaceAngularFluxes(Write("no_upwind", without_upwind), {"right"}, mpi_comm),
    std::logic_error);
  // Interior tags have no upwind values.
  EXPECT_NO_THROW(DistributeSurfaceAngularFluxes(
    Write("no_upwind_interior", without_upwind), {"plane"}, mpi_comm));
}

TEST_F(SurfaceResponseDataTest, RejectsFileOfAnotherSet)
{
  // File 2 of the five was left over from another set written to the same prefix.
  auto spec = ForwardSpec();
  spec.stray_file = std::make_pair(size_t{2}, std::uint64_t{8});
  EXPECT_THROW(DistributeSurfaceAngularFluxes(Write("stray", spec), all_tags, mpi_comm),
               std::logic_error);
}

TEST_F(SurfaceResponseDataTest, RejectsQuadratureDifferingAcrossFiles)
{
  // File 3 of the five has other quadrature weights, so it was not written with the rest.
  auto spec = ForwardSpec();
  spec.odd_quadrature_file = 3;
  EXPECT_THROW(DistributeSurfaceAngularFluxes(Write("odd", spec), all_tags, mpi_comm),
               std::logic_error);
}

TEST_F(SurfaceResponseDataTest, RejectsDifferentGroupRanges)
{
  // Same group counts per groupset, but the adjoint groupsets start one group later.
  auto spec = AdjointSpec();
  spec.first_group_shift = 1;
  const auto forward =
    DistributeSurfaceAngularFluxes(Write("forward", ForwardSpec()), all_tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(Write("adjoint", spec), all_tags, mpi_comm);
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"right"}, mpi_comm),
               std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, RejectsMissingFace)
{
  const auto forward_base = Write("forward", ForwardSpec());
  auto spec = AdjointSpec();
  spec.tags[1].faces_per_file = {{}, {0, 1, 2, 3, 5, 6}, {}};
  const auto adjoint_base = Write("adjoint", spec);

  const auto forward = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(adjoint_base, all_tags, mpi_comm);
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"plane"}, mpi_comm), std::logic_error);
  // The boundary is complete.
  EXPECT_NO_THROW(ComputeSurfaceResponses(forward, adjoint, {"right"}, mpi_comm));
}

TEST_F(SurfaceResponseDataTest, RejectsDuplicateFace)
{
  auto spec = ForwardSpec();
  spec.tags[0].faces_per_file = {{0, 1, 2, 3, 4}, {2}, {5}, {}, {6}};
  const auto forward_base = Write("forward", spec);
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm), std::logic_error);
}

TEST_F(SurfaceResponseDataTest, RejectsMissingTagAndFiles)
{
  const auto forward_base = Write("forward", ForwardSpec());
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, {"absent"}, mpi_comm),
               std::logic_error);
  EXPECT_THROW(
    DistributeSurfaceAngularFluxes((directory_ / "absent_p").string(), all_tags, mpi_comm),
    std::logic_error);
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, {}, mpi_comm), std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, RejectsIncompatibleDataSets)
{
  const auto forward_base = Write("forward", ForwardSpec());
  auto spec = AdjointSpec();
  spec.centroid_spacing = 2.0e-9;
  const auto adjoint_base = Write("adjoint", spec);

  const auto forward = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(adjoint_base, all_tags, mpi_comm);
  // Keys of different spacings cannot be matched.
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm), std::invalid_argument);
  // Forward data in the adjoint role.
  EXPECT_THROW(ComputeSurfaceResponses(forward, forward, {}, mpi_comm), std::invalid_argument);
  // Unknown surface.
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"left"}, mpi_comm),
               std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, MatchesDirectionsRelativeToQuadratureSpacing)
{
  const auto forward =
    DistributeSurfaceAngularFluxes(Write("forward", ForwardSpec()), all_tags, mpi_comm);
  const auto adjoint =
    DistributeSurfaceAngularFluxes(Write("adjoint", AdjointSpec()), all_tags, mpi_comm);
  // The nearest two fixture directions are 1.2 apart, so directions match within 1.2e-6.
  const auto Shift = [](DistributedSurfaceData data, const std::uint64_t d, const double delta)
  {
    for (auto& quadrature : data.groupsets)
      quadrature.omegas[3 * d] += delta;
    return data;
  };
  const double interior = ExpectedResponse(true);
  const double boundary = ExpectedResponse(false);

  // Within the tolerance: every direction, stored in reverse order, still finds its match.
  auto close = adjoint;
  for (auto& quadrature : close.groupsets)
    for (auto& component : quadrature.omegas)
      component += 1.0e-9;
  const auto surface_responses = ComputeSurfaceResponses(forward, close, {}, mpi_comm);
  ASSERT_EQ(surface_responses.size(), 2U);
  EXPECT_NEAR(surface_responses[0].outgoing + surface_responses[0].incoming,
              interior,
              1.0e-13 * std::abs(interior));
  EXPECT_NEAR(surface_responses[1].Total(), boundary, 1.0e-13 * std::abs(boundary));

  // Beyond the tolerance.
  EXPECT_THROW(ComputeSurfaceResponses(forward, Shift(adjoint, 0, 1.2e-4), {}, mpi_comm),
               std::invalid_argument);

  // Duplicate adjoint directions.
  auto duplicate = adjoint;
  for (auto& quadrature : duplicate.groupsets)
    std::copy_n(quadrature.omegas.begin(), 3, quadrature.omegas.begin() + 3);
  EXPECT_THROW(ComputeSurfaceResponses(forward, duplicate, {}, mpi_comm), std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, MatchesCloselySpacedDirections)
{
  auto forward =
    DistributeSurfaceAngularFluxes(Write("forward", ForwardSpec()), all_tags, mpi_comm);
  auto adjoint =
    DistributeSurfaceAngularFluxes(Write("adjoint", AdjointSpec()), all_tags, mpi_comm);
  // Place fixture direction 1 1e-12 from direction 0 in both sets. The adjoint stores directions in
  // reverse order, so there they are adjoint directions 2 and 3. Only the direction map reads the
  // vectors, so the closed-form response is unchanged.
  const auto Place = [](DistributedSurfaceData& data, const std::uint64_t d, const std::uint64_t to)
  {
    for (auto& quadrature : data.groupsets)
    {
      std::copy_n(quadrature.omegas.begin() + static_cast<std::ptrdiff_t>(3 * to),
                  3,
                  quadrature.omegas.begin() + static_cast<std::ptrdiff_t>(3 * d));
      quadrature.omegas[3 * d + 2] += 1.0e-12;
    }
  };
  Place(forward, 1, 0);
  Place(adjoint, 2, 3);

  const auto surface_responses = ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm);
  ASSERT_EQ(surface_responses.size(), 2U);
  const double interior = ExpectedResponse(true);
  const double boundary = ExpectedResponse(false);
  EXPECT_NEAR(surface_responses[0].outgoing + surface_responses[0].incoming,
              interior,
              1.0e-13 * std::abs(interior));
  EXPECT_NEAR(surface_responses[1].Total(), boundary, 1.0e-13 * std::abs(boundary));

  // At this spacing, a 1e-15 difference is 1e-3 of the spacing and no longer a match.
  const std::size_t shifted_direction = 2;
  for (auto& quadrature : adjoint.groupsets)
    quadrature.omegas[3 * shifted_direction] += 1.0e-15;
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm), std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, RejectsMismatchedWeights)
{
  auto adjoint_spec = AdjointSpec();
  adjoint_spec.weight_scale = 2.0;
  const auto forward =
    DistributeSurfaceAngularFluxes(Write("forward", ForwardSpec()), all_tags, mpi_comm);
  const auto adjoint =
    DistributeSurfaceAngularFluxes(Write("adjoint", adjoint_spec), all_tags, mpi_comm);
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm), std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, ReaderReturnsOneFileOfTheSet)
{
  const auto forward_base = Write("forward", ForwardSpec());
  // Every rank reads the same file; the reader is not collective.
  const auto surf_fluxes =
    DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, all_tags);
  // Groupset-major: one entry per groupset and surface.
  ASSERT_EQ(surf_fluxes.size(), groups_per_groupset.size() * all_tags.size());
  for (size_t gs = 0; gs < groups_per_groupset.size(); ++gs)
    for (size_t s = 0; s < all_tags.size(); ++s)
    {
      const auto& surf_flux = surf_fluxes[gs * all_tags.size() + s];
      EXPECT_EQ(surf_flux.groupset_id, static_cast<int>(gs));
      EXPECT_EQ(surf_flux.surface_name, all_tags[s]);
      EXPECT_FALSE(surf_flux.adjoint);
      EXPECT_EQ(surf_flux.num_partitions, 5U);
      EXPECT_EQ(surf_flux.num_directions, directions.size());
      EXPECT_EQ(surf_flux.num_groups, groups_per_groupset[gs]);
      const auto num_faces = ForwardSpec().tags[s].faces_per_file[0].size();
      EXPECT_EQ(surf_flux.mapping->cell_ids.size(), num_faces);
      EXPECT_EQ(surf_flux.mapping->cell_map.size(), num_faces);
      EXPECT_EQ(surf_flux.mapping->normals.size(), 3 * num_faces);
      EXPECT_EQ(surf_flux.mapping->mass_matrix.size(), 4 * num_faces);
      // The face geometry is shared by the groupsets of a tag.
      EXPECT_EQ(surf_flux.mapping, surf_fluxes[s].mapping);
      EXPECT_EQ(surf_flux.data.omega.size(), 3 * directions.size());
      EXPECT_EQ(surf_flux.data.weights.size(), directions.size());
      EXPECT_EQ(surf_flux.data.psi.size(),
                2 * num_faces * directions.size() * surf_flux.num_groups);
    }
  EXPECT_THROW(DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, {}),
               std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, RejectsInconsistentDatasetSizes)
{
  auto spec = ForwardSpec();
  spec.truncate_flux = true;
  const auto forward_base = Write("truncated", spec);
  EXPECT_THROW(DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, all_tags),
               std::logic_error);
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm), std::logic_error);
}

TEST_F(SurfaceResponseDataTest, RejectsNonFiniteData)
{
  auto spec = ForwardSpec();
  spec.nan_flux = true;
  const auto forward_base = Write("nan", spec);
  EXPECT_THROW(DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, all_tags),
               std::logic_error);
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm), std::logic_error);
}

TEST_F(SurfaceResponseDataTest, RejectsPartitionOutsideItsSet)
{
  // Every file records a set of two files, so files 2 to 4 cannot belong to it.
  auto spec = ForwardSpec();
  spec.recorded_partitions = 2;
  const auto forward_base = Write("misnamed", spec);
  EXPECT_NO_THROW(DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 1, all_tags));
  EXPECT_THROW(DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 2, all_tags),
               std::logic_error);
}

TEST_F(SurfaceResponseDataTest, ResolvesInteriorSurfaceNames)
{
  const auto forward_base = Write("forward", ForwardSpec());

  // The reader reads both sides of a named interior surface, `_u` first, per groupset.
  const auto surf_fluxes =
    DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, {"right", "plane"});
  ASSERT_EQ(surf_fluxes.size(), 3 * groups_per_groupset.size());
  EXPECT_EQ(surf_fluxes[0].surface_name, "right");
  EXPECT_EQ(surf_fluxes[1].surface_name, "plane_u");
  EXPECT_EQ(surf_fluxes[2].surface_name, "plane_d");

  // Distributing by name gives the same data as the stored tags.
  const auto by_name = DistributeSurfaceAngularFluxes(forward_base, {"plane", "right"}, mpi_comm);
  const auto by_tag = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm);
  ASSERT_EQ(by_name.surfaces.size(), by_tag.surfaces.size());
  for (const auto& tag : all_tags)
  {
    ASSERT_EQ(by_name.surfaces.count(tag), 1U);
    EXPECT_EQ(by_name.surfaces.at(tag).index, by_tag.surfaces.at(tag).index);
    EXPECT_EQ(by_name.surfaces.at(tag).coords, by_tag.surfaces.at(tag).coords);
  }

  // Unknown names, names with only one stored side, and names that resolve to the same tag twice
  // are rejected.
  EXPECT_THROW(DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, {"absent"}),
               std::logic_error);
  EXPECT_THROW(
    DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(forward_base, 0, {"plane", "plane_u"}),
    std::invalid_argument);
  auto one_side = ForwardSpec();
  one_side.tags.erase(one_side.tags.begin() + 1);
  EXPECT_THROW(DistributeSurfaceAngularFluxes(Write("one_side", one_side), {"plane"}, mpi_comm),
               std::logic_error);
}

TEST_F(SurfaceResponseDataTest, InteriorPartsMatchClosedForm)
{
  const auto forward =
    DistributeSurfaceAngularFluxes(Write("forward", ForwardSpec()), {"plane"}, mpi_comm);
  const auto adjoint =
    DistributeSurfaceAngularFluxes(Write("adjoint", AdjointSpec()), {"plane"}, mpi_comm);
  const auto responses = ComputeSurfaceResponses(forward, adjoint, {"plane"}, mpi_comm);
  ASSERT_EQ(responses.size(), 1U);
  // Directions crossing along +x pair the forward `_u` side with the adjoint `_d` side, and those
  // crossing back the forward `_d` side with the adjoint `_u` side.
  const auto [outgoing, incoming] = ExpectedInteriorParts();
  EXPECT_NEAR(responses.front().outgoing, outgoing, 1.0e-13 * std::abs(outgoing));
  EXPECT_NEAR(responses.front().incoming, incoming, 1.0e-13 * std::abs(incoming));
  EXPECT_DOUBLE_EQ(responses.front().Total(), outgoing + incoming);
}

TEST_F(SurfaceResponseDataTest, RejectsLoneInteriorSide)
{
  // A surface response needs both sides of an interior surface, in either data set.
  const auto forward_base = Write("forward", ForwardSpec());
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, {"plane_u"}, mpi_comm),
               std::invalid_argument);
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, {"plane_d", "right"}, mpi_comm),
               std::invalid_argument);
  EXPECT_NO_THROW(DistributeSurfaceAngularFluxes(forward_base, {"plane_u", "plane_d"}, mpi_comm));
}

TEST_F(SurfaceResponseDataTest, SideTagsOrientInteriorSurfaces)
{
  const auto forward_base = Write("forward", ForwardSpec());
  const auto adjoint_base = Write("adjoint", AdjointSpec());
  const auto forward = DistributeSurfaceAngularFluxes(forward_base, all_tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(adjoint_base, all_tags, mpi_comm);
  const auto interior = ComputeSurfaceResponses(forward, adjoint, {"plane"}, mpi_comm);
  ASSERT_EQ(interior.size(), 1U);
  const double outgoing = interior.front().outgoing;
  const double incoming = interior.front().incoming;

  // With both sides in the forward data, a side tag selects the interior surface along that side's
  // normal: `_u` is the default, and `_d` reverses every cosine, so its outgoing directions are
  // those leaving the `_d` cell and both terms change sign.
  const auto plane_u = ComputeSurfaceResponses(forward, adjoint, {"plane_u"}, mpi_comm);
  ASSERT_EQ(plane_u.size(), 1U);
  EXPECT_EQ(plane_u.front().name, "plane_u");
  EXPECT_EQ(plane_u.front().kind, SurfaceKind::interior);
  EXPECT_EQ(plane_u.front().outgoing, outgoing);
  EXPECT_EQ(plane_u.front().incoming, incoming);
  const auto plane_d = ComputeSurfaceResponses(forward, adjoint, {"plane_d"}, mpi_comm);
  ASSERT_EQ(plane_d.size(), 1U);
  EXPECT_EQ(plane_d.front().name, "plane_d");
  EXPECT_EQ(plane_d.front().kind, SurfaceKind::interior);
  EXPECT_EQ(plane_d.front().outgoing, -incoming);
  EXPECT_EQ(plane_d.front().incoming, -outgoing);

  // A surface may be requested once under any of its names, and only interior surfaces take a
  // side tag.
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"plane", "plane_d"}, mpi_comm),
               std::invalid_argument);
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"plane_u", "plane_d"}, mpi_comm),
               std::invalid_argument);
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"right_d"}, mpi_comm),
               std::invalid_argument);
}

TEST_F(SurfaceResponseDataTest, BoundaryTagsWithSideSuffixesStayBoundaries)
{
  // Boundaries named like interior tags: `outer_d` alone, and `side_u` and `side_d`, which look
  // like both sides of an interior surface. Each holds the data of the `right` boundary.
  const std::vector<std::string> tags = {"outer_d", "side_u", "side_d"};
  const auto RenameBoundary = [&tags](SetSpec spec)
  {
    const auto right = spec.tags.back();
    spec.tags.clear();
    for (const auto& name : tags)
    {
      auto tag = right;
      tag.name = name;
      spec.tags.push_back(tag);
    }
    return spec;
  };
  const auto forward_base = Write("forward", RenameBoundary(ForwardSpec()));
  const auto adjoint_base = Write("adjoint", RenameBoundary(AdjointSpec()));
  const auto forward = DistributeSurfaceAngularFluxes(forward_base, tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(adjoint_base, tags, mpi_comm);

  // The recorded kinds, not the names, classify the tags.
  const auto surface_responses = ComputeSurfaceResponses(forward, adjoint, {}, mpi_comm);
  ASSERT_EQ(surface_responses.size(), 3U);
  const double boundary = ExpectedResponse(false);
  for (const auto& surface_response : surface_responses)
  {
    EXPECT_EQ(surface_response.kind, SurfaceKind::boundary) << surface_response.name;
    EXPECT_NEAR(surface_response.Total(), boundary, 1.0e-13 * std::abs(boundary));
  }
  EXPECT_NO_THROW(ComputeSurfaceResponses(forward, adjoint, {"outer_d"}, mpi_comm));

  // A bare name reads two tags only when both are interior sides.
  EXPECT_THROW(DistributeSurfaceAngularFluxes(forward_base, {"side"}, mpi_comm), std::logic_error);
}

TEST_F(SurfaceResponseDataTest, RejectsAdjointTagsOfAnotherKind)
{
  // The adjoint set records the interior tags as boundaries.
  auto spec = AdjointSpec();
  for (auto& tag : spec.tags)
    if (tag.side != 'b')
      tag.recorded_interior = false;
  const auto forward =
    DistributeSurfaceAngularFluxes(Write("forward", ForwardSpec()), all_tags, mpi_comm);
  const auto adjoint = DistributeSurfaceAngularFluxes(Write("adjoint", spec), all_tags, mpi_comm);
  EXPECT_THROW(ComputeSurfaceResponses(forward, adjoint, {"plane"}, mpi_comm),
               std::invalid_argument);
  EXPECT_NO_THROW(ComputeSurfaceResponses(forward, adjoint, {"right"}, mpi_comm));
}

TEST(SurfaceFaceOwnerTest, IsInRangeAndDeterministic)
{
  for (int num_ranks : {1, 2, 3, 7})
    for (std::int64_t k = -5; k < 5; ++k)
    {
      const DiscreteOrdinatesProblemIO::FaceCentroidKey key = {k, 3 * k, -k};
      const auto owner = SurfaceFaceOwner(key, num_ranks);
      EXPECT_GE(owner, 0);
      EXPECT_LT(owner, num_ranks);
      EXPECT_EQ(owner, SurfaceFaceOwner(key, num_ranks));
    }
}
