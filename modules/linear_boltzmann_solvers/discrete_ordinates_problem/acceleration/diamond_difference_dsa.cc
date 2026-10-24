// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/diamond_difference_dsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <array>
#include <cmath>

namespace opensn
{

struct DiamondDifferenceDSA::GroupSystem
{
  Mat matrix = nullptr;
  Vec rhs = nullptr;
  Vec solution = nullptr;
  KSP solver = nullptr;

  GroupSystem() = default;
  GroupSystem(const GroupSystem&) = delete;
  GroupSystem& operator=(const GroupSystem&) = delete;

  ~GroupSystem()
  {
    KSPDestroy(&solver);
    VecDestroy(&solution);
    VecDestroy(&rhs);
    MatDestroy(&matrix);
  }
};

DiamondDifferenceDSA::~DiamondDifferenceDSA()
{
  systems_.clear();
  if (communicator_ != MPI_COMM_NULL)
    MPI_Comm_free(&communicator_);
}

DiamondDifferenceDSA::DiamondDifferenceDSA(DiscreteOrdinatesProblem& problem) : problem_(problem)
{
  const auto& groupset = problem_.GetGroupsets().front();
  OpenSnInvalidArgumentIf(groupset.angle_agg->GetNumDelayedAngularDOFs().second != 0,
                          "Diamond-difference DSA does not support lagged interface fluxes.");
  edges_ = PieceWiseLinearContinuous::New(problem_.GetGrid());
  unknowns_.AddUnknown(UnknownType::SCALAR);
  const auto local_size = static_cast<PetscInt>(edges_->GetNumLocalDOFs(unknowns_));
  const auto global_size = static_cast<PetscInt>(edges_->GetNumGlobalDOFs(unknowns_));
  const auto ghost_ids = edges_->GetGhostDOFIndices(unknowns_);
  const std::vector<PetscInt> ghosts(ghost_ids.begin(), ghost_ids.end());
  std::vector<int64_t> diagonal_nnz, off_diagonal_nnz;
  edges_->BuildSparsityPattern(diagonal_nnz, off_diagonal_nnz, unknowns_);

  double first_half_moment = 0.0;
  double second_half_moment = 0.0;
  for (std::size_t n = 0; n < groupset.quadrature->GetNumAngles(); ++n)
  {
    const double mu = groupset.quadrature->GetOmega(n).z;
    if (mu > 0.0)
    {
      const double weight = groupset.quadrature->GetWeight(n);
      first_half_moment += weight * mu;
      second_half_moment += weight * mu * mu;
    }
  }
  OpenSnInvalidArgumentIf(second_half_moment <= 0.0,
                          "Diamond-difference DSA requires a nonzero half-range second moment.");
  const double alpha = first_half_moment / (3.0 * second_half_moment);

  OpenSnMPICall(MPI_Comm_dup(opensn::mpi_comm, &communicator_));
  systems_.reserve(groupset.GetNumGroups());
  for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
  {
    auto system = std::make_unique<GroupSystem>();
    system->matrix = CreateSquareMatrix(communicator_, local_size, global_size);
    InitMatrixSparsity(system->matrix, diagonal_nnz, off_diagonal_nnz);
    OpenSnPETScCall(VecCreateGhost(communicator_,
                                   local_size,
                                   global_size,
                                   static_cast<PetscInt>(ghosts.size()),
                                   ghosts.data(),
                                   &system->rhs));
    OpenSnPETScCall(VecDuplicate(system->rhs, &system->solution));

    for (const auto& cell : problem_.GetGrid()->GetLocalCells())
    {
      const auto& mapping = edges_->GetCellMapping(*cell);
      const auto& nodes = mapping.GetNodeLocations();
      const double h = (nodes[1] - nodes[0]).Norm();
      const auto& xs = *problem_.GetBlockID2XSMap().at(cell->block_id);
      const double sigma_t = xs.GetSigmaTotal()[g];
      const double sigma_r = sigma_t - xs.GetTransferMatrix(0).GetValueIJ(g, g);
      const double stiffness = 1.0 / (3.0 * sigma_t * h);
      const double mass = sigma_r * h / 4.0;
      const std::array<PetscInt, 2> rows{
        static_cast<PetscInt>(edges_->MapDOF(*cell, 0, unknowns_, 0, 0)),
        static_cast<PetscInt>(edges_->MapDOF(*cell, 1, unknowns_, 0, 0))};
      const std::array<double, 4> values{
        stiffness + mass, -stiffness + mass, -stiffness + mass, stiffness + mass};
      OpenSnPETScCall(
        MatSetValues(system->matrix, 2, rows.data(), 2, rows.data(), values.data(), ADD_VALUES));
      for (std::size_t f = 0; f < cell->faces.size(); ++f)
        if (not cell->faces[f].has_neighbor)
        {
          const auto node = mapping.MapFaceNode(f, 0);
          OpenSnPETScCall(MatSetValue(system->matrix, rows[node], rows[node], alpha, ADD_VALUES));
        }
    }
    OpenSnPETScCall(MatAssemblyBegin(system->matrix, MAT_FINAL_ASSEMBLY));
    OpenSnPETScCall(MatAssemblyEnd(system->matrix, MAT_FINAL_ASSEMBLY));
    OpenSnPETScCall(KSPCreate(communicator_, &system->solver));
    OpenSnPETScCall(KSPSetOperators(system->solver, system->matrix, system->matrix));
    const std::string prefix = problem_.GetName() + "_dd_dsa_g" + std::to_string(g) + "_";
    OpenSnPETScCall(KSPSetOptionsPrefix(system->solver, prefix.c_str()));
    PC pc = nullptr;
    OpenSnPETScCall(KSPGetPC(system->solver, &pc));
    OpenSnPETScCall(KSPSetType(system->solver, KSPCG));
    OpenSnPETScCall(PCSetType(pc, PCGAMG));
    OpenSnPETScCall(KSPSetTolerances(system->solver, 1.0e-12, PETSC_DEFAULT, PETSC_DEFAULT, 1000));
    OpenSnPETScCall(KSPSetErrorIfNotConverged(system->solver, PETSC_TRUE));
    OpenSnPETScCall(KSPSetFromOptions(system->solver));
    systems_.push_back(std::move(system));
  }
}

void
DiamondDifferenceDSA::Apply()
{
  CALI_CXX_MARK_SCOPE("DiamondDifferenceDSA::Apply");
  auto& phi = problem_.GetPhiNewLocal();
  const auto& old_phi = problem_.GetPhiOldLocal();
  std::vector<double> delta(phi.size());
  for (std::size_t i = 0; i < phi.size(); ++i)
    delta[i] = phi[i] - old_phi[i];
  std::vector<double> previous_correction(edges_->GetNumLocalAndGhostDOFs(unknowns_), 0.0);
  std::vector<double> correction;
  const auto& transport = problem_.GetSpatialDiscretization();
  const auto& phi_unknowns = problem_.GetUnknownManager();

  for (unsigned int g = 0; g < systems_.size(); ++g)
  {
    auto& system = *systems_[g];
    OpenSnPETScCall(VecSet(system.rhs, 0.0));
    for (const auto& cell : problem_.GetGrid()->GetLocalCells())
    {
      const auto& nodes = edges_->GetCellMapping(*cell).GetNodeLocations();
      const double h = (nodes[1] - nodes[0]).Norm();
      const auto& scatter = problem_.GetBlockID2XSMap().at(cell->block_id)->GetTransferMatrix(0);
      double source_sum = 0.0;
      for (unsigned int i = 0; i < 2; ++i)
      {
        source_sum +=
          scatter.GetValueIJ(g, g) * delta[transport.MapDOFLocal(*cell, i, phi_unknowns, 0, g)];
        if (g > 0)
          source_sum += scatter.GetValueIJ(g, g - 1) *
                        (delta[transport.MapDOFLocal(*cell, i, phi_unknowns, 0, g - 1)] +
                         previous_correction[edges_->MapDOFLocal(*cell, i, unknowns_, 0, 0)]);
      }
      for (unsigned int i = 0; i < 2; ++i)
      {
        const auto row = static_cast<PetscInt>(edges_->MapDOF(*cell, i, unknowns_, 0, 0));
        OpenSnPETScCall(VecSetValue(system.rhs, row, h * source_sum / 4.0, ADD_VALUES));
      }
    }
    OpenSnPETScCall(VecAssemblyBegin(system.rhs));
    OpenSnPETScCall(VecAssemblyEnd(system.rhs));
    OpenSnPETScCall(VecSet(system.solution, 0.0));
    OpenSnPETScCall(KSPSolve(system.solver, system.rhs, system.solution));
    CommunicateGhostEntries(system.solution);
    edges_->LocalizePETScVectorWithGhosts(system.solution, correction, unknowns_);
    for (const auto& cell : problem_.GetGrid()->GetLocalCells())
      for (unsigned int i = 0; i < 2; ++i)
        phi[transport.MapDOFLocal(*cell, i, phi_unknowns, 0, g)] +=
          correction[edges_->MapDOFLocal(*cell, i, unknowns_, 0, 0)];
    previous_correction.swap(correction);
  }
}

} // namespace opensn
