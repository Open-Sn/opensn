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

// "side" and "face_x2" are the same operator on orthogonal hexahedra (non-uniform spacing).
TEST(DiffusionMIPSidePenalty, EqualsFaceX2OnOrthogonalHexes)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  auto grid = BuildOrthogonalMesh({{0.0, 0.5, 0.75, 1.5}, {0.0, 0.2, 1.0}, {0.0, 0.3, 0.4}});
  grid->SetUniformBlockID(0);
  grid->SetOrthogonalBoundaries();
  for (const auto& bc : std::vector<BoundaryCondition>{{BCType::DIRICHLET, {0, 0, 0}},
                                                       {BCType::ROBIN, {0.25, 0.5, 0}}})
  {
    PetscInt n_side = 0, n_face = 0;
    const auto A_side = MIPMatrix(grid, "side", bc, n_side);
    const auto A_face = MIPMatrix(grid, "face_x2", bc, n_face);
    ASSERT_EQ(n_side, n_face);
    double diff = 0.0, scale = 0.0;
    for (size_t k = 0; k < A_side.size(); ++k)
    {
      diff = std::max(diff, std::fabs(A_side[k] - A_face[k]));
      scale = std::max(scale, std::fabs(A_face[k]));
    }
    EXPECT_LT(diff, 1.0e-12 * scale);
  }
}

// Severely distorted quadrilaterals (Shestakov's random mesh, parameter a = 0.05, 8 x 8 cells of
// the unit square): the polygon penalty of the other conventions (length V/A_f, no factor 2) leaves
// the operator indefinite, so DSA diverges; "side" keeps it positive definite with either boundary
// form.
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
    PetscInt n = 0;
    EXPECT_FALSE(IsPositiveDefinite(MIPMatrix(grid, "cell", bc, n), n));
    EXPECT_TRUE(IsPositiveDefinite(MIPMatrix(grid, "side", bc, n), n));
  }
}
