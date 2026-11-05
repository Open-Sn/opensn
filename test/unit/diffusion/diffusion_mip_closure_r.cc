// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gmock/gmock.h"
#include "test/unit/common/mesh_builders.h"

#include "modules/diffusion/boundary_closure.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/acceleration.h"

#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
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

/// Whether a symmetric dense matrix (row-major n x n) is positive definite: its Cholesky
/// factorization exists with positive pivots.
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

/// MIP matrix of the unit cube (num_cells^3 hexes) with every boundary set to `bc`, one group with
/// D = 1 / (3 sigma_t) and sigma_a = 1e-4 sigma_t, and the closure quadrature `quad`.
std::vector<double>
CubeMatrix(unsigned int num_cells,
           double sigma_t,
           const BoundaryCondition& bc,
           const AngularQuadrature& quad,
           PetscInt& n)
{
  auto grid = BuildBoxMesh(1.0, num_cells, 0.0);
  grid->SetUniformBlockID(0);
  grid->SetOrthogonalBoundaries();
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  std::map<uint64_t, BoundaryCondition> bcs;
  for (const auto& [bid, name] : grid->GetBoundaryIDMap())
    bcs[bid] = bc;
  std::map<unsigned int, Multigroup_D_and_sigR> xs;
  xs.insert(std::make_pair(0, Multigroup_D_and_sigR{{1.0 / (3.0 * sigma_t)}, {1.0e-4 * sigma_t}}));
  std::vector<UnitCellMatrices> unit_cell_matrices(grid->GetLocalCellCount());
  for (const auto& cell : grid->GetLocalCells())
    unit_cell_matrices[cell->local_id] = ComputeUnitCellIntegrals(sdm, *cell);

  DiffusionMIPSolver solver(
    "ClosureR", sdm, sdm.UNITARY_UNKNOWN_MANAGER, bcs, xs, unit_cell_matrices, false, false);
  solver.options.perform_symmetry_check = true;
  solver.options.solver_policy = "direct";
  solver.SetBoundaryClosureQuadrature(quad.GetOmegas(), quad.GetWeights());
  solver.Initialize();
  std::vector<double> q(sdm.GetNumLocalDOFs(sdm.UNITARY_UNKNOWN_MANAGER), 0.0);
  solver.AssembleAand_b(q);
  return DenseCopy(solver.GetSystemMatrix(), n);
}

} // namespace

// beta_f of the 8-direction GLC 2x4 set: every direction has |Omega . n| = 1/sqrt(3) on the axis
// faces, and on the (1,1,1) diagonal two directions have 1 and six have 1/3 (mean 1/2). A fine set
// approaches the exact 1/4.
TEST(DiffusionClosureR, QuadratureHalfRangeMoment)
{
  GLCProductQuadrature3DXYZ s2(2, 4, 0);
  const double inv_sqrt3 = 1.0 / std::sqrt(3.0);
  for (const Vector3 n : {Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)})
    EXPECT_NEAR(
      QuadratureHalfRangeMoment(s2.GetOmegas(), s2.GetWeights(), n), 0.5 * inv_sqrt3, 1.0e-12);
  EXPECT_NEAR(QuadratureHalfRangeMoment(
                s2.GetOmegas(), s2.GetWeights(), Vector3(inv_sqrt3, inv_sqrt3, inv_sqrt3)),
              0.25,
              1.0e-12);

  GLCProductQuadrature3DXYZ fine(32, 64, 0);
  for (const Vector3 n :
       {Vector3(1, 0, 0), Vector3(0, 0, 1), Vector3(inv_sqrt3, inv_sqrt3, inv_sqrt3)})
    EXPECT_NEAR(QuadratureHalfRangeMoment(fine.GetOmegas(), fine.GetWeights(), n), 0.25, 1.0e-3);
}

// The blend: 1 on thin cells, (s0 / sigma_t h)^2 on thick ones, and where the coercivity root binds
// it solves (1 - t)^2 / (1 + t) = 6 beta sigma_t h / (C_grad n).
TEST(DiffusionClosureR, Theta)
{
  EXPECT_DOUBLE_EQ(ClosureRTheta(0.25, 0.5, 1.0, 1, 1.5), 1.0);
  EXPECT_DOUBLE_EQ(ClosureRTheta(0.25, 1.5, 1.0, 1, 1.5), 1.0);
  EXPECT_NEAR(ClosureRTheta(0.25, 1000.0, 1.0, 1, 1.5), std::pow(1.5 / 1000.0, 2), 1.0e-15);

  // PWQD-like trace constant on a corner cell: the smooth rule gives 0.25, the root 0.455.
  const double beta = 0.25, sth = 3.0, c_grad = 7.36;
  const unsigned int n = 3;
  const double theta = ClosureRTheta(beta, sth, c_grad, n, 1.5);
  const double r = 6.0 * beta * sth / (c_grad * n);
  EXPECT_GT(theta, 0.25);
  EXPECT_NEAR((1.0 - theta) * (1.0 - theta) / (1.0 + theta), r, 1.0e-12);

  // A cell without diffusion (void-like D = 0): sigma_t h is infinite, theta 0.
  EXPECT_DOUBLE_EQ(ClosureRTheta(0.25, std::numeric_limits<double>::infinity(), 1.0, 1, 1.5), 0.0);
}

// Trace constant of PWLD: 1 on the faces of a cube (the prototype's value, snproto
// MIPDiffusion._trace_constant), exactly 3 on tetrahedra (linear functions: (g.n)^2 A_f / (|g|^2 V)
// with h_f = 3V/A_f).
TEST(DiffusionClosureR, TraceConstant)
{
  {
    auto grid = BuildBoxMesh(2.0, 2, 0.0);
    std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto K = ComputeUnitCellIntegrals(*sdm, *cell).intV_gradshapeI_gradshapeJ;
      for (unsigned int f = 0; f < cell->faces.size(); ++f)
        EXPECT_NEAR(FaceTraceConstant(sdm->GetCellMapping(*cell), *cell, f, K), 1.0, 1.0e-10);
    }
  }
  {
    const auto mfile =
      std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / "gmsh_all_tets.vtu";
    auto grid = BuildMeshFromFile(mfile);
    std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
    size_t checked = 0;
    for (const auto& cell : grid->GetLocalCells())
    {
      if (checked++ == 50)
        break;
      const auto K = ComputeUnitCellIntegrals(*sdm, *cell).intV_gradshapeI_gradshapeJ;
      for (unsigned int f = 0; f < cell->faces.size(); ++f)
        EXPECT_NEAR(FaceTraceConstant(sdm->GetCellMapping(*cell), *cell, f, K), 3.0, 1.0e-9);
    }
  }
}

// In RZ both norms of the trace constant carry the weight r, so the constant does not depend on the
// distance to the axis: far from it the weighted constant is the Cartesian one. (With an unweighted
// face matrix against the weighted stiffness it would scale like 1/r.)
TEST(DiffusionClosureR, TraceConstantRZ)
{
  const double r0 = 1000.0;
  auto grid = BuildOrthogonalMesh({{r0, r0 + 0.5, r0 + 1.0}, {0.0, 0.3, 1.0}});
  std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& mapping = sdm->GetCellMapping(*cell);
    const auto K = ComputeUnitCellIntegrals(*sdm, *cell).intV_gradshapeI_gradshapeJ;
    const auto K_rz = ComputeUnitCellIntegrals(*sdm, *cell, CoordinateSystemType::CYLINDRICAL)
                        .intV_gradshapeI_gradshapeJ;
    for (unsigned int f = 0; f < cell->faces.size(); ++f)
    {
      const double c_xy = FaceTraceConstant(mapping, *cell, f, K);
      const double c_rz =
        FaceTraceConstant(mapping, *cell, f, K_rz, CoordinateSystemType::CYLINDRICAL);
      EXPECT_NEAR(c_rz, c_xy, 2.0e-3 * c_xy);
    }
  }
}

// The assembled operator with every boundary face closed by closure R is symmetric (symmetry check
// at assembly) and positive definite on thin, intermediate and thick cells.
TEST(DiffusionClosureR, OperatorIsSPD)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  // The check itself rejects an indefinite matrix (eigenvalues 3 and -1).
  EXPECT_FALSE(IsPositiveDefinite({1.0, 2.0, 2.0, 1.0}, 2));

  GLCProductQuadrature3DXYZ quad(2, 4, 0);
  const unsigned int num_cells = 3;
  for (const double sigma_t_h : {0.01, 1.0, 1000.0})
  {
    SCOPED_TRACE(sigma_t_h);
    PetscInt n = 0;
    const auto A =
      CubeMatrix(num_cells, sigma_t_h * num_cells, {BCType::CLOSURE_R, {0, 0, 0}}, quad, n);
    EXPECT_TRUE(IsPositiveDefinite(A, n));
  }
}

// On thin cells theta = 1 and closure R is the natural Marshak form with the quadrature's
// coefficient 2 beta_f: the same matrix as ROBIN {2 beta_f, 1, 0}. On thick cells it differs.
TEST(DiffusionClosureR, EqualsRobinOnThinCells)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  GLCProductQuadrature3DXYZ quad(2, 4, 0);
  const double two_beta = 1.0 / std::sqrt(3.0); // 2 beta_f on the axis faces of GLC 2x4
  const unsigned int num_cells = 3;
  for (const double sigma_t_h : {0.01, 1000.0})
  {
    PetscInt n = 0, m = 0;
    const auto R =
      CubeMatrix(num_cells, sigma_t_h * num_cells, {BCType::CLOSURE_R, {0, 0, 0}}, quad, n);
    const auto B =
      CubeMatrix(num_cells, sigma_t_h * num_cells, {BCType::ROBIN, {two_beta, 1.0, 0.0}}, quad, m);
    ASSERT_EQ(n, m);
    double diff = 0.0, scale = 0.0;
    for (size_t k = 0; k < R.size(); ++k)
    {
      diff = std::max(diff, std::fabs(R[k] - B[k]));
      scale = std::max(scale, std::fabs(B[k]));
    }
    if (sigma_t_h < 1.0)
      EXPECT_LT(diff, 1.0e-12 * scale);
    else
      EXPECT_GT(diff, 1.0e-3 * scale);
  }
}
