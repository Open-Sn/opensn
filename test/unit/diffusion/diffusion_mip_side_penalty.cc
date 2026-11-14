// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gmock/gmock.h"
#include "test/unit/common/mesh_builders.h"

#include "modules/diffusion/boundary_closure.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"

#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"

#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

/// Dense copy of an assembled PETSc matrix (single rank), row-major.
std::vector<double>
DenseCopy(const Mat& A, PetscInt& n)
{
  PetscInt m = 0;
  MatGetSize(A, &n, &m);
  std::vector<double> dense(static_cast<size_t>(n * n), 0.0);
  for (PetscInt i = 0; i < n; ++i)
  {
    PetscInt ncols = 0;
    const PetscInt* cols = nullptr;
    const PetscScalar* vals = nullptr;
    MatGetRow(A, i, &ncols, &cols, &vals);
    for (PetscInt k = 0; k < ncols; ++k)
      dense[i * n + cols[k]] = vals[k];
    MatRestoreRow(A, i, &ncols, &cols, &vals);
  }
  return dense;
}

/// Whether a symmetric dense matrix (row-major n x n) is positive definite (Cholesky pivots > 0).
bool
IsPositiveDefinite(std::vector<double> a, PetscInt n)
{
  for (PetscInt k = 0; k < n; ++k)
  {
    double pivot = a[k * n + k];
    for (PetscInt p = 0; p < k; ++p)
      pivot -= a[k * n + p] * a[k * n + p];
    if (not(pivot > 0.0))
      return false;
    const double l_kk = std::sqrt(pivot);
    a[k * n + k] = l_kk;
    for (PetscInt i = k + 1; i < n; ++i)
    {
      double v = a[i * n + k];
      for (PetscInt p = 0; p < k; ++p)
        v -= a[i * n + p] * a[k * n + p];
      a[i * n + k] = v / l_kk;
    }
  }
  return true;
}

/// Dense MIP matrix of `grid` (one block, every boundary set to `bc`) for an h_perp convention, on
/// thin cells (D = 100, sigma_a = 1e-3) where the penalty alone has to keep the operator definite.
std::vector<double>
MIPMatrix(const std::shared_ptr<Mesh>& grid,
          const std::string& convention,
          const BoundaryCondition& bc,
          PetscInt& n)
{
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;
  std::map<uint64_t, BoundaryCondition> bcs;
  for (const auto& [bid, name] : grid->GetBoundaryIDMap())
    bcs[bid] = bc;
  std::map<unsigned int, Multigroup_D_and_sigR> xs;
  xs.insert(std::make_pair(0, Multigroup_D_and_sigR{{100.0}, {1.0e-3}}));
  std::vector<UnitCellMatrices> unit_cell_matrices(grid->GetLocalCellCount());
  for (const auto& cell : grid->GetLocalCells())
    unit_cell_matrices[cell->local_id] = ComputeUnitCellIntegrals(sdm, *cell);
  DiffusionMIPSolver solver(
    "SidePenalty", sdm, sdm.UNITARY_UNKNOWN_MANAGER, bcs, xs, unit_cell_matrices, false, false);
  solver.options.h_perp_convention = convention;
  solver.options.solver_policy = "direct";
  solver.Initialize();
  std::vector<double> q(sdm.GetNumLocalDOFs(sdm.UNITARY_UNKNOWN_MANAGER), 0.0);
  solver.AssembleAand_b(q);
  return DenseCopy(solver.GetSystemMatrix(), n);
}

} // namespace

// On rectangles and boxes the side triangles / tetrahedra of a face reach the cell centre, half the
// cell width away, so the side length 2 min_k h_k is the cell width normal to the face, V/A_f.
TEST(DiffusionMIPSidePenalty, SideLengthIsWidthOnOrthogonalCells)
{
  const std::vector<double> x = {0.0, 0.5, 0.75, 1.5}, y = {0.0, 0.2, 1.0}, z = {0.0, 0.3, 0.4};
  for (const auto& node_sets : std::vector<std::vector<std::vector<double>>>{{x, y}, {x, y, z}})
  {
    auto grid = BuildOrthogonalMesh(node_sets);
    for (const auto& cell : grid->GetLocalCells())
      for (unsigned int f = 0; f < cell->faces.size(); ++f)
        EXPECT_NEAR(FaceSideLength(*grid, *cell, f), cell->volume / cell->faces[f].area, 1.0e-12);
  }
}

// A cell centroid on the line (2D) or plane (3D) of a face flattens that face's side simplices;
// FaceSideLength rejects it instead of returning a zero length (an infinite penalty).
TEST(DiffusionMIPSidePenalty, DegenerateSideSimplexThrows)
{
  const std::vector<double> x = {0.0, 0.5, 1.0}, y = {0.0, 0.2, 1.0}, z = {0.0, 0.3, 0.4};
  for (const auto& node_sets : std::vector<std::vector<std::vector<double>>>{{x, y}, {x, y, z}})
  {
    auto grid = BuildOrthogonalMesh(node_sets);
    const auto& cell = *grid->GetLocalCells().front();
    Cell flat = cell;
    flat.centroid = flat.faces[0].centroid;
    EXPECT_NO_THROW(FaceSideLength(*grid, cell, 0));
    EXPECT_THROW(FaceSideLength(*grid, flat, 0), std::logic_error);
  }
}

// "face_side" = min(face length, 2 x side length): on rectangles and boxes the cap is inactive and
// the length is the cell width normal to the face, V/A_f, with the factor 2 on polygons.
TEST(DiffusionMIPSidePenalty, FaceSideIsFaceLengthOnOrthogonalCells)
{
  const std::vector<double> x = {0.0, 0.5, 0.75, 1.5}, y = {0.0, 0.2, 1.0}, z = {0.0, 0.3, 0.4};
  for (const auto& node_sets : std::vector<std::vector<std::vector<double>>>{{x, y}, {x, y, z}})
  {
    auto grid = BuildOrthogonalMesh(node_sets);
    grid->SetUniformBlockID(0);
    std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
    const auto& sdm = *sdm_ptr;
    std::map<uint64_t, BoundaryCondition> bcs;
    std::map<unsigned int, Multigroup_D_and_sigR> xs;
    xs.insert(std::make_pair(0, Multigroup_D_and_sigR{{1.0}, {1.0}}));
    std::vector<UnitCellMatrices> unit_cell_matrices(grid->GetLocalCellCount());
    DiffusionMIPSolver solver(
      "FaceSide", sdm, sdm.UNITARY_UNKNOWN_MANAGER, bcs, xs, unit_cell_matrices, false, false);
    solver.options.h_perp_convention = "face_side";
    for (const auto& cell : grid->GetLocalCells())
      for (unsigned int f = 0; f < cell->faces.size(); ++f)
        EXPECT_NEAR(solver.HPerpendicular(*cell, f), cell->volume / cell->faces[f].area, 1.0e-12);
    EXPECT_DOUBLE_EQ(solver.PolygonPenaltyScale(), 2.0);
  }
}

// Severely distorted quadrilaterals (Shestakov's random mesh, parameter a = 0.05, 8 x 8 cells of
// the unit square): OpenSn's original polygon penalty ("cell": length V/A_f, no factor 2) leaves
// the operator indefinite, so DSA diverges; "face_side" keeps it positive definite with either
// boundary form, because the side-length cap engages on the faces with a skinny side triangle.
TEST(DiffusionMIPSidePenalty, PositiveDefiniteOnDistortedQuads)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto mfile =
    std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / "shestakov_quad_a0.05_n8.vtu";
  auto grid = BuildMeshFromFile(mfile);
  grid->SetUniformBlockID(0);
  grid->SetOrthogonalBoundaries();
  for (const auto& bc : std::vector<BoundaryCondition>{{BCType::DIRICHLET, {0, 0, 0}},
                                                       {BCType::ROBIN, {0.25, 0.5, 0}}})
  {
    // Build each matrix before the call: the order in which function arguments are evaluated is
    // unspecified, so passing n alongside MIPMatrix(..., n) could pass n before it is set.
    PetscInt n_cell = 0, n_face_side = 0;
    const auto A_cell = MIPMatrix(grid, "cell", bc, n_cell);
    const auto A_face_side = MIPMatrix(grid, "face_side", bc, n_face_side);
    ASSERT_GT(n_cell, 0);
    ASSERT_EQ(n_face_side, n_cell);
    EXPECT_FALSE(IsPositiveDefinite(A_cell, n_cell));
    EXPECT_TRUE(IsPositiveDefinite(A_face_side, n_face_side));
  }
  // The cap is active somewhere: twice the side length is below the face length on some face.
  double min_ratio = 1.0;
  for (const auto& cell : grid->GetLocalCells())
    for (unsigned int f = 0; f < cell->faces.size(); ++f)
      min_ratio =
        std::min(min_ratio, 2.0 * FaceSideLength(*grid, *cell, f) / FaceOrthogonalLength(*cell, f));
  EXPECT_LT(min_ratio, 0.1);
}

// On every tetrahedron the cell centroid is a quarter of the height above each face, so the side
// length is half of 3V/A_f and "face_side" is the face length 3V/A_f itself (the cap is inactive);
// the same on cubes, where both lengths are V/A_f.
TEST(DiffusionMIPSidePenalty, FaceSideIsFaceLengthOnTetsAndCubes)
{
  const auto mfile =
    std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / "gmsh_all_tets.vtu";
  for (const auto& grid : {BuildMeshFromFile(mfile), BuildBoxMesh(1.0, 3, 0.0)})
  {
    grid->SetUniformBlockID(0);
    std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
    const auto& sdm = *sdm_ptr;
    std::map<uint64_t, BoundaryCondition> bcs;
    std::map<unsigned int, Multigroup_D_and_sigR> xs;
    xs.insert(std::make_pair(0, Multigroup_D_and_sigR{{1.0}, {1.0}}));
    std::vector<UnitCellMatrices> unit_cell_matrices(grid->GetLocalCellCount());
    DiffusionMIPSolver solver(
      "FaceSide", sdm, sdm.UNITARY_UNKNOWN_MANAGER, bcs, xs, unit_cell_matrices, false, false);
    solver.options.h_perp_convention = "face_side";
    for (const auto& cell : grid->GetLocalCells())
      for (unsigned int f = 0; f < cell->faces.size(); ++f)
      {
        const double h_face = FaceOrthogonalLength(*cell, f);
        EXPECT_NEAR(solver.HPerpendicular(*cell, f), h_face, 1.0e-10 * h_face);
        if (cell->faces.size() == 4)
          EXPECT_NEAR(FaceSideLength(*grid, *cell, f), 0.5 * h_face, 1.0e-10 * h_face);
      }
  }
}

// Groupset defaults: closure R and the "face_side" length for WGDSA and TGDSA.
TEST(DiffusionMIPSidePenalty, GroupsetDefaults)
{
  auto params = LBSGroupset::GetInputParameters();
  EXPECT_EQ(params.GetParamValue<std::string>("wgdsa_boundary_closure"), "closure_r");
  EXPECT_EQ(params.GetParamValue<std::string>("tgdsa_boundary_closure"), "closure_r");
  EXPECT_EQ(params.GetParamValue<std::string>("wgdsa_h_perp_convention"), "face_side");
  EXPECT_EQ(params.GetParamValue<std::string>("tgdsa_h_perp_convention"), "face_side");
}
