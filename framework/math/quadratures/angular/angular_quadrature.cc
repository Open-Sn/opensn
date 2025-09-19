// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <cmath>

#include <petscmat.h>
#include <petscksp.h>

namespace opensn
{

void
AngularQuadrature::MakeHarmonicIndices()
{
  m_to_ell_em_map_.clear();

  if (dimension_ == 1)
    for (auto ell = 0; ell <= scattering_order_; ++ell)
      m_to_ell_em_map_.emplace_back(ell, 0);
  else if (dimension_ == 2)
    for (auto ell = 0; ell <= scattering_order_; ++ell)
      for (auto m = -ell; m <= ell; m += 2)
        m_to_ell_em_map_.emplace_back(ell, m);
  else if (dimension_ == 3)
    for (auto ell = 0; ell <= scattering_order_; ++ell)
      for (auto m = -ell; m <= ell; ++m)
        m_to_ell_em_map_.emplace_back(ell, m);
}

std::vector<std::vector<double>>
AngularQuadrature::InvertMatrix(const std::vector<std::vector<double>>& matrix)
{
  const size_t n = matrix.size();

  // Check if matrix is square
  if (n == 0 || matrix[0].size() != n)
  {
    throw std::runtime_error("Matrix must be square for inversion");
  }

  PetscErrorCode ierr;
  Mat A, B, X;
  KSP ksp;

  // Create PETSc matrix A from input matrix
  ierr = MatCreate(PETSC_COMM_SELF, &A);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetSizes(A, n, n, n, n);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetType(A, MATSEQDENSE);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetUp(A);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Fill matrix A with input data
  for (PetscInt i = 0; i < static_cast<PetscInt>(n); ++i)
  {
    for (PetscInt j = 0; j < static_cast<PetscInt>(n); ++j)
    {
      ierr = MatSetValue(A, i, j, matrix[i][j], INSERT_VALUES);
      CHKERRABORT(PETSC_COMM_SELF, ierr);
    }
  }
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Create identity matrix B (right-hand side)
  ierr = MatCreate(PETSC_COMM_SELF, &B);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetSizes(B, n, n, n, n);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetType(B, MATSEQDENSE);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetUp(B);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Fill B with identity matrix
  for (PetscInt i = 0; i < static_cast<PetscInt>(n); ++i)
  {
    ierr = MatSetValue(B, i, i, 1.0, INSERT_VALUES);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
  }
  ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Create solution matrix X
  ierr = MatDuplicate(B, MAT_DO_NOT_COPY_VALUES, &X);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Set up linear solver
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPSetOperators(ksp, A, A);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Use direct solver for better accuracy
  PC pc;
  ierr = KSPGetPC(ksp, &pc);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = PCSetType(pc, PCLU);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Solve AX = B (where B is identity, so X will be A^-1)
  ierr = KSPMatSolve(ksp, B, X);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Check convergence
  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  if (reason < 0)
  {
    // Clean up before throwing
    MatDestroy(&A);
    MatDestroy(&B);
    MatDestroy(&X);
    KSPDestroy(&ksp);
    throw std::runtime_error("PETSc linear solver failed to converge during matrix inversion");
  }

  // Extract solution back to std::vector format
  std::vector<std::vector<double>> inverse(n, std::vector<double>(n));

  // Get array from PETSc dense matrix
  PetscScalar* array;
  ierr = MatDenseGetArray(X, &array);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Copy data (PETSc uses column-major storage for dense matrices)
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      inverse[i][j] = PetscRealPart(array[j * n + i]); // Transpose due to column-major
    }
  }

  ierr = MatDenseRestoreArray(X, &array);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  // Clean up PETSc objects
  ierr = MatDestroy(&A);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatDestroy(&B);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatDestroy(&X);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  log.Log0Verbose1() << "Matrix inversion completed successfully using PETSc LU solver";

  return inverse;
}

void
AngularQuadrature::BuildDiscreteToMomentOperatorStandard()
{
  d2m_op_.clear();

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  for (const auto& ell_em : m_to_ell_em_map_)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (auto n = 0; n < num_angles; ++n)
    {
      const auto& cur_angle = abscissae[n];
      double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
      double w = weights[n];
      cur_mom.push_back(value * w);
    }

    d2m_op_.push_back(cur_mom);
  }

  // Verbose printout
  std::stringstream outs;
  outs << "\nQuadrature d2m operator (Standard Method):\n";
  for (auto n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (auto m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << d2m_op_[m][n]
           << " ";
    }
    outs << "\n";
  }

  log.Log0Verbose1() << outs.str();
}

void
AngularQuadrature::BuildMomentToDiscreteOperatorStandard()
{
  m2d_op_.clear();

  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  const auto normalization = std::accumulate(weights.begin(), weights.end(), 0.0);

  for (const auto& ell_em : m_to_ell_em_map_)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (auto n = 0; n < num_angles; ++n)
    {
      const auto& cur_angle = abscissae[n];
      double value = ((2.0 * ell_em.ell + 1.0) / normalization) *
                     Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
      cur_mom.push_back(value);
    }

    m2d_op_.push_back(cur_mom);
  }

  // Verbose printout
  std::stringstream outs;

  outs << "\nQuadrature m2d operator (Standard Method):\n";
  for (auto n = 0; n < num_angles; ++n)
  {
    outs << std::setw(5) << n;
    for (auto m = 0; m < num_moms; ++m)
    {
      outs << std::setw(15) << std::left << std::fixed << std::setprecision(10) << m2d_op_[m][n]
           << " ";
    }
    outs << "\n";
  }

  log.Log0Verbose1() << outs.str();
}

void
AngularQuadrature::BuildDiscreteToMomentOperator()
{
  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  switch (construction_method_)
  {
    case OperatorConstructionMethod::Standard:
    {
      BuildDiscreteToMomentOperatorStandard();
      break;
    }

    case OperatorConstructionMethod::GalerkinMethodOne:
    {
      log.Log0Verbose1() << "Building D2M operator by inverting M2D operator using PETSc";

      // Build M2D first if not already built
      if (m2d_op_.empty())
      {
        BuildMomentToDiscreteOperatorStandard();
      }

      // Check dimensions
      if (num_angles != num_moms)
      {
        throw std::runtime_error("Cannot invert M2D operator: number of directions (" +
                                 std::to_string(num_angles) + ") != number of moments (" +
                                 std::to_string(num_moms) + ")");
      }

      try
      {
        d2m_op_ = InvertMatrix(m2d_op_);
        log.Log0Verbose1() << "D2M operator successfully computed as inverse of M2D using PETSc";
      }
      catch (const std::exception& e)
      {
        log.LogAllError() << "Failed to invert M2D operator with PETSc: " << e.what();
        log.LogAllError() << "Falling back to Standard Method";
        BuildDiscreteToMomentOperatorStandard();
      }
      break;
    }

    case OperatorConstructionMethod::GalerkinMethodTwo:
    {
      BuildDiscreteToMomentOperatorStandard();
      break;
    }
  }
}

void
AngularQuadrature::BuildMomentToDiscreteOperator()
{
  const size_t num_angles = abscissae.size();
  const size_t num_moms = m_to_ell_em_map_.size();

  switch (construction_method_)
  {
    case OperatorConstructionMethod::Standard:
    {
      BuildMomentToDiscreteOperatorStandard();
      break;
    }

    case OperatorConstructionMethod::GalerkinMethodOne:
    {
      BuildMomentToDiscreteOperatorStandard();
      break;
    }

    case OperatorConstructionMethod::GalerkinMethodTwo:
    {
      log.Log0Verbose1() << "Building M2D operator by inverting D2M operator using PETSc";

      // Build D2M first if not already built
      if (d2m_op_.empty())
      {
        BuildDiscreteToMomentOperatorStandard();
      }

      // Check dimensions
      if (num_angles != num_moms)
      {
        throw std::runtime_error("Cannot invert D2M operator: number of directions (" +
                                 std::to_string(num_angles) + ") != number of moments (" +
                                 std::to_string(num_moms) + ")");
      }

      try
      {
        m2d_op_ = InvertMatrix(d2m_op_);
        log.Log0Verbose1() << "M2D operator successfully computed as inverse of D2M using PETSc";
      }
      catch (const std::exception& e)
      {
        log.LogAllError() << "Failed to invert D2M operator with PETSc: " << e.what();
        log.LogAllError() << "Falling back to Standard Method";
        BuildMomentToDiscreteOperatorStandard();
      }
      break;
    }
  }
}

std::vector<std::vector<double>> const&
AngularQuadrature::GetDiscreteToMomentOperator() const
{
  return d2m_op_;
}

std::vector<std::vector<double>> const&
AngularQuadrature::GetMomentToDiscreteOperator() const
{
  return m2d_op_;
}

const std::vector<AngularQuadrature::HarmonicIndices>&
AngularQuadrature::GetMomentToHarmonicsIndexMap() const
{
  return m_to_ell_em_map_;
}

} // namespace opensn
