// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/cbc_spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds_adams_adams_hawkins.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_set/aah_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/source_functions/source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/logging/log.h"
#include "framework/logging/log_exceptions.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <iomanip>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, DiscreteOrdinatesSolver);

DiscreteOrdinatesSolver::DiscreteOrdinatesSolver(const std::string& name) : LBSSolver(name)
{
}

InputParameters
DiscreteOrdinatesSolver::GetInputParameters()
{
  InputParameters params = LBSSolver::GetInputParameters();

  params.SetClassName("DiscreteOrdinatesSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "LBSDiscreteOrdinatesSolver");

  params.AddOptionalParameterArray(
    "directions_sweep_order_to_print",
    std::vector<size_t>(),
    "List of direction id's for which sweep ordering info is to be printed.");

  params.AddOptionalParameter(
    "sweep_type", "AAH", "The sweep type to use for sweep operatorations.");

  params.ConstrainParameterRange("sweep_type", AllowableRangeList::New({"AAH", "CBC"}));

  return params;
}

DiscreteOrdinatesSolver::DiscreteOrdinatesSolver(const InputParameters& params)
  : LBSSolver(params),
    verbose_sweep_angles_(params.GetParamVectorValue<size_t>("directions_sweep_order_to_print")),
    sweep_type_(params.GetParamValue<std::string>("sweep_type"))
{
}

DiscreteOrdinatesSolver::~DiscreteOrdinatesSolver()
{
  CALI_CXX_MARK_FUNCTION;

  for (auto& groupset : groupsets_)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    // Reset sweep orderings
    if (groupset.angle_agg_ != nullptr)
      groupset.angle_agg_->angle_set_groups.clear();
  }
}

std::pair<size_t, size_t>
DiscreteOrdinatesSolver::GetNumPhiIterativeUnknowns()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::GetNumPhiIterativeUnknowns");
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  size_t num_local_psi_dofs = 0;
  size_t num_globl_psi_dofs = 0;
  for (auto& groupset : groupsets_)
  {
    const auto num_delayed_psi_info = groupset.angle_agg_->GetNumDelayedAngularDOFs();
    num_local_psi_dofs += num_delayed_psi_info.first;
    num_globl_psi_dofs += num_delayed_psi_info.second;
  }

  const size_t num_local_dofs = num_local_phi_dofs + num_local_psi_dofs;
  const size_t num_globl_dofs = num_globl_phi_dofs + num_globl_psi_dofs;

  return {num_local_dofs, num_globl_dofs};
}

void
DiscreteOrdinatesSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::Initialize");

  LBSSolver::Initialize();

  // Initialize source func
  using namespace std::placeholders;
  auto src_function = std::make_shared<SourceFunction>(*this);
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);

  // Initialize groupsets for sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }
  InitializeSolverSchemes();
}

void
DiscreteOrdinatesSolver::InitializeWGSSolvers()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::InitializeWGSSolvers");

  wgs_solvers_.clear(); // this is required
  for (auto& groupset : groupsets_)
  {
    std::shared_ptr<SweepChunk> sweep_chunk = SetSweepChunk(groupset);
    auto sweep_wgs_context_ptr = std::make_shared<SweepWGSContext>(
      *this,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options_.verbose_inner_iterations,
      sweep_chunk);

    if (groupset.iterative_method_ == IterativeMethod::CLASSIC_RICHARDSON)
      wgs_solvers_.push_back(std::make_shared<ClassicRichardson>(sweep_wgs_context_ptr));
    else
      wgs_solvers_.push_back(std::make_shared<WGSLinearSolver>(sweep_wgs_context_ptr));
  }
}

void
DiscreteOrdinatesSolver::ScalePhiVector(PhiSTLOption which_phi, double value)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::ScalePhiVector");

  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  Scale(*y_ptr, value);

  for (auto& groupset : groupsets_)
  {
    switch (which_phi)
    {
      case PhiSTLOption::PHI_NEW:
      {
        auto psi = groupset.angle_agg_->GetNewDelayedAngularDOFsAsSTLVector();
        Scale(psi, value);
        groupset.angle_agg_->SetNewDelayedAngularDOFsFromSTLVector(psi);
        break;
      }
      case PhiSTLOption::PHI_OLD:
      {
        auto psi = groupset.angle_agg_->GetOldDelayedAngularDOFsAsSTLVector();
        Scale(psi, value);
        groupset.angle_agg_->SetOldDelayedAngularDOFsFromSTLVector(psi);
        break;
      }
    }
  }
}

void
DiscreteOrdinatesSolver::SetGSPETScVecFromPrimarySTLvector(const LBSGroupset& groupset,
                                                           Vec x,
                                                           PhiSTLOption which_phi)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::SetGSPETScVecFromPrimarySTLvector");

  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x, &x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          x_ref[index] = (*y_ptr)[mapping + g]; // Offset on purpose
        }                                       // for g
      }                                         // for moment
    }                                           // for dof
  }                                             // for cell

  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      groupset.angle_agg_->AppendNewDelayedAngularDOFsToArray(index, x_ref);
      break;
    case PhiSTLOption::PHI_OLD:
      groupset.angle_agg_->AppendOldDelayedAngularDOFsToArray(index, x_ref);
      break;
  }

  VecRestoreArray(x, &x_ref);
}

void
DiscreteOrdinatesSolver::SetPrimarySTLvectorFromGSPETScVec(const LBSGroupset& groupset,
                                                           Vec x,
                                                           PhiSTLOption which_phi)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::SetPrimarySTLvectorFromGSPETScVec");

  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromGSPETScVec");
  }

  const double* x_ref;
  VecGetArrayRead(x, &x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          (*y_ptr)[mapping + g] = x_ref[index];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell

  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      groupset.angle_agg_->SetNewDelayedAngularDOFsFromArray(index, x_ref);
      break;
    case PhiSTLOption::PHI_OLD:
      groupset.angle_agg_->SetOldDelayedAngularDOFsFromArray(index, x_ref);
  }

  VecRestoreArrayRead(x, &x_ref);
}

void
DiscreteOrdinatesSolver::GSScopedCopyPrimarySTLvectors(const LBSGroupset& groupset,
                                                       PhiSTLOption from_which_phi,
                                                       PhiSTLOption to_which_phi)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::GSScopedCopyPrimarySTLvectors");

  std::vector<double>* y_ptr;
  switch (to_which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  std::vector<double>* x_src_ptr;
  switch (from_which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      x_src_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      x_src_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  int gsi = groupset.groups_.front().id_;
  size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          (*y_ptr)[mapping + g] = (*x_src_ptr)[mapping + g];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell

  if (from_which_phi == PhiSTLOption::PHI_NEW and to_which_phi == PhiSTLOption::PHI_OLD)
    groupset.angle_agg_->SetDelayedPsiOld2New();
  if (from_which_phi == PhiSTLOption::PHI_OLD and to_which_phi == PhiSTLOption::PHI_NEW)
    groupset.angle_agg_->SetDelayedPsiNew2Old();
}

void
DiscreteOrdinatesSolver::SetMultiGSPETScVecFromPrimarySTLvector(
  const std::vector<int>& groupset_ids, Vec x, PhiSTLOption which_phi)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::SetMultiGSPETScVecFromPrimarySTLvector");

  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetMultiGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x, &x_ref);

  int64_t index = -1;
  for (int gs_id : groupset_ids)
  {
    auto& groupset = groupsets_.at(gs_id);

    int gsi = groupset.groups_.front().id_;
    int gsf = groupset.groups_.back().id_;
    int gss = gsf - gsi + 1;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      auto& transport_view = cell_transport_views_[cell.local_id_];

      for (int i = 0; i < cell.vertex_ids_.size(); i++)
      {
        for (int m = 0; m < num_moments_; m++)
        {
          size_t mapping = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; g++)
          {
            index++;
            x_ref[index] = (*y_ptr)[mapping + g]; // Offset on purpose
          }                                       // for g
        }                                         // for moment
      }                                           // for dof
    }                                             // for cell

    switch (which_phi)
    {
      case PhiSTLOption::PHI_NEW:
        groupset.angle_agg_->AppendNewDelayedAngularDOFsToArray(index, x_ref);
        break;
      case PhiSTLOption::PHI_OLD:
        groupset.angle_agg_->AppendOldDelayedAngularDOFsToArray(index, x_ref);
        break;
    }
  } // for groupset id

  VecRestoreArray(x, &x_ref);
}

void
DiscreteOrdinatesSolver::SetPrimarySTLvectorFromMultiGSPETScVecFrom(
  const std::vector<int>& groupset_ids, Vec x, PhiSTLOption which_phi)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::SetPrimarySTLvectorFromMultiGSPETScVecFrom");

  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromMultiGSPETScVecFrom");
  }

  const double* x_ref;
  VecGetArrayRead(x, &x_ref);

  int64_t index = -1;
  for (int gs_id : groupset_ids)
  {
    auto& groupset = groupsets_.at(gs_id);

    int gsi = groupset.groups_.front().id_;
    int gsf = groupset.groups_.back().id_;
    int gss = gsf - gsi + 1;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      auto& transport_view = cell_transport_views_[cell.local_id_];

      for (int i = 0; i < cell.vertex_ids_.size(); i++)
      {
        for (int m = 0; m < num_moments_; m++)
        {
          size_t mapping = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; g++)
          {
            index++;
            (*y_ptr)[mapping + g] = x_ref[index];
          } // for g
        }   // for moment
      }     // for dof
    }       // for cell

    switch (which_phi)
    {
      case PhiSTLOption::PHI_NEW:
        groupset.angle_agg_->SetNewDelayedAngularDOFsFromArray(index, x_ref);
        break;
      case PhiSTLOption::PHI_OLD:
        groupset.angle_agg_->SetOldDelayedAngularDOFsFromArray(index, x_ref);
    }
  } // for groupset id

  VecRestoreArrayRead(x, &x_ref);
}

void
DiscreteOrdinatesSolver::ReorientAdjointSolution()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::ReorientAdjointSolution");

  for (const auto& groupset : groupsets_)
  {
    int gs = groupset.id_;

    // Moment map for flux moments
    const auto& moment_map = groupset.quadrature_->GetMomentToHarmonicsIndexMap();

    // Angular flux info
    auto& psi = psi_new_local_[gs];
    const auto& uk_man = groupset.psi_uk_man_;

    // Build reversed angle mapping
    std::map<int, int> reversed_angle_map;
    if (options_.save_angular_flux)
    {
      const auto& omegas = groupset.quadrature_->omegas_;
      const auto num_gs_angles = omegas.size();

      // Go through angles until all are paired
      std::set<size_t> visited;
      for (int idir = 0; idir < num_gs_angles; ++idir)
      {
        // Skip if already encountered
        if (visited.count(idir) > 0)
          continue;

        bool found = true;
        for (int jdir = 0; jdir < num_gs_angles; ++jdir)
        {
          // Angles are opposite if their sum is zero
          const auto sum = grid_ptr_->Dimension() == 1
                             ? Vector3(0.0, 0.0, omegas[idir].z + omegas[jdir].z)
                             : omegas[idir] + omegas[jdir];
          const bool opposite = sum.NormSquare() < 1.0e-8;

          // Add opposites to mapping
          if (opposite)
          {
            found = true;
            reversed_angle_map[idir] = jdir;

            visited.insert(idir);
            visited.insert(jdir);
          }
        } // for angle n

        OpenSnLogicalErrorIf(not found,
                             "Opposing angle for " + omegas[idir].PrintStr() + " in groupset " +
                               std::to_string(gs) + " not found.");

      } // for angle m
    }   // if saving angular flux

    const auto num_gs_groups = groupset.groups_.size();
    const auto gsg_i = groupset.groups_.front().id_;
    const auto gsg_f = groupset.groups_.back().id_;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& transport_view = cell_transport_views_[cell.local_id_];
      for (int i = 0; i < transport_view.NumNodes(); ++i)
      {
        // Reorient flux moments
        //
        // Because flux moments are integrated angular fluxes, the
        // angular flux and spherical harmonics must be evaluated at
        // opposite angles in the quadrature integration. Taking advantage
        // of the even/odd nature of the spherical harmonics, i.e.
        // Y_{\ell,m}(-\Omega) = (-1)^\ell Y_{\ell,m}(\Omega), the flux
        // moments must be multiplied by (-1)^\ell.
        for (int imom = 0; imom < num_moments_; ++imom)
        {
          const auto& ell = moment_map[imom].ell;
          const auto dof_map = transport_view.MapDOF(i, imom, 0);

          for (int g = gsg_i; g <= gsg_f; ++g)
          {
            phi_new_local_[dof_map + g] *= std::pow(-1.0, ell);
            phi_old_local_[dof_map + g] *= std::pow(-1.0, ell);
          } // for group g
        }   // for moment m

        // Reorient angular flux
        if (options_.save_angular_flux)
        {
          for (const auto& [idir, jdir] : reversed_angle_map)
          {
            const auto dof_map =
              std::make_pair(discretization_->MapDOFLocal(cell, i, uk_man, idir, 0),
                             discretization_->MapDOFLocal(cell, i, uk_man, jdir, 0));

            for (int gsg = 0; gsg < num_gs_groups; ++gsg)
              std::swap(psi[dof_map.first + gsg], psi[dof_map.second + gsg]);
          }
        }
      } // for node i
    }   // for cell

  } // for groupset
}

void
DiscreteOrdinatesSolver::ZeroOutflowBalanceVars(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::ZeroOutflowBalanceVars");

  for (const auto& cell : grid_ptr_->local_cells)
    for (int f = 0; f < cell.faces_.size(); ++f)
      for (auto& group : groupset.groups_)
        cell_transport_views_[cell.local_id_].ZeroOutflow(f, group.id_);
}

void
DiscreteOrdinatesSolver::ComputeBalance()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::ComputeBalance");

  opensn::mpi_comm.barrier();

  // Get material source
  // This is done using the SetSource routine because it allows a lot of flexibility.
  auto mat_src = phi_old_local_;
  mat_src.assign(mat_src.size(), 0.0);
  for (auto& groupset : groupsets_)
  {
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    active_set_source_function_(groupset,
                                q_moments_local_,
                                phi_old_local_,
                                APPLY_FIXED_SOURCES | APPLY_AGS_FISSION_SOURCES |
                                  APPLY_WGS_FISSION_SOURCES);
    LBSSolver::GSScopedCopyPrimarySTLvectors(groupset, q_moments_local_, mat_src);
  }

  // Compute absorption, material-source and in-flow
  double local_out_flow = 0.0;
  double local_in_flow = 0.0;
  double local_absorption = 0.0;
  double local_production = 0.0;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = discretization_->GetCellMapping(cell);
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& fe_intgrl_values = unit_cell_matrices_[cell.local_id_];
    const size_t num_nodes = transport_view.NumNodes();
    const auto& IntV_shapeI = fe_intgrl_values.intV_shapeI;
    const auto& IntS_shapeI = fe_intgrl_values.intS_shapeI;

    // Inflow: This is essentially an integration over all faces, all angles, and all groups. For
    // non-reflective boundaries, only the cosines that are negative are added to the inflow
    // integral. For reflective boundaries, it is expected that, upon convergence, inflow = outflow
    // (within numerical tolerances set by the user).
    for (int f = 0; f < cell.faces_.size(); ++f)
    {
      const auto& face = cell.faces_[f];

      if (not face.has_neighbor_) // Boundary face
      {
        const auto& bndry = sweep_boundaries_[face.neighbor_id_];

        if (bndry->IsReflecting())
        {
          for (int g = 0; g < num_groups_; ++g)
            local_in_flow += transport_view.GetOutflow(f, g);
        }
        else
        {
          for (const auto& groupset : groupsets_)
          {
            for (int n = 0; n < groupset.quadrature_->omegas_.size(); ++n)
            {
              const auto& omega = groupset.quadrature_->omegas_[n];
              const double wt = groupset.quadrature_->weights_[n];
              const double mu = omega.Dot(face.normal_);

              if (mu < 0.0)
              {
                for (int fi = 0; fi < face.vertex_ids_.size(); ++fi)
                {
                  const int i = cell_mapping.MapFaceNode(f, fi);
                  const auto& IntFi_shapeI = IntS_shapeI[f][i];

                  for (const auto& group : groupset.groups_)
                  {
                    const int g = group.id_;
                    const double psi = *bndry->PsiIncoming(cell.local_id_, f, fi, n, g, 0);
                    local_in_flow -= mu * wt * psi * IntFi_shapeI;
                  } // for group
                }   // for fi
              }     // if mu < 0
            }       // for n
          }         // for groupset
        }           // if reflecting boundary
      }             // if boundary
    }               // for f

    // Outflow: The group-wise outflow was determined during a solve so we just accumulate it here.
    for (int f = 0; f < cell.faces_.size(); ++f)
    {
      const auto& face = cell.faces_[f];
      for (int g = 0; g < num_groups_; ++g)
        local_out_flow += transport_view.GetOutflow(f, g);
    }

    // Absorption and sources
    const auto& xs = transport_view.XS();
    const auto& sigma_a = xs.SigmaAbsorption();
    for (int i = 0; i < num_nodes; ++i)
    {
      for (int g = 0; g < num_groups_; ++g)
      {
        size_t imap = transport_view.MapDOF(i, 0, g);
        double phi_0g = phi_old_local_[imap];
        double q_0g = mat_src[imap];

        local_absorption += sigma_a[g] * phi_0g * IntV_shapeI[i];
        local_production += q_0g * IntV_shapeI[i];
      } // for g
    }   // for i
  }     // for cell

  // Compute local balance
  double local_balance = local_production + local_in_flow - local_absorption - local_out_flow;
  double local_gain = local_production + local_in_flow;
  std::vector<double> local_balance_table = {
    local_absorption, local_production, local_in_flow, local_out_flow, local_balance, local_gain};
  size_t table_size = local_balance_table.size();

  // Compute global balance
  std::vector<double> globl_balance_table(table_size, 0.0);
  mpi_comm.all_reduce(
    local_balance_table.data(), table_size, globl_balance_table.data(), mpi::op::sum<double>());
  double globl_absorption = globl_balance_table.at(0);
  double globl_production = globl_balance_table.at(1);
  double globl_in_flow = globl_balance_table.at(2);
  double globl_out_flow = globl_balance_table.at(3);
  double globl_balance = globl_balance_table.at(4);
  double globl_gain = globl_balance_table.at(5);

  log.Log() << "Balance table:\n"
            << std::setprecision(6) << std::scientific
            << " Absorption rate             = " << globl_absorption << "\n"
            << " Production rate             = " << globl_production << "\n"
            << " In-flow rate                = " << globl_in_flow << "\n"
            << " Out-flow rate               = " << globl_out_flow << "\n"
            << " Gain (In-flow + Production) = " << globl_gain << "\n"
            << " Balance (Gain - Loss)       = " << globl_balance << "\n"
            << " Balance/Gain, in %          = " << globl_balance / globl_gain * 100. << "\n";

  opensn::mpi_comm.barrier();
}

std::vector<double>
DiscreteOrdinatesSolver::ComputeLeakage(const unsigned int groupset_id,
                                        const uint64_t boundary_id) const
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::ComputeLeakage");

  // Perform checks
  OpenSnInvalidArgumentIf(groupset_id < 0 or groupset_id >= groupsets_.size(),
                          "Invalid groupset id.");
  OpenSnLogicalErrorIf(not options_.save_angular_flux,
                       "The option `save_angular_flux` must be set to `true` in order "
                       "to compute outgoing currents.");

  const auto& sdm = *discretization_;
  const auto& groupset = groupsets_.at(groupset_id);
  const auto& psi_uk_man = groupset.psi_uk_man_;
  const auto& quadrature = groupset.quadrature_;

  const auto num_gs_angles = quadrature->omegas_.size();
  const auto num_gs_groups = groupset.groups_.size();

  const auto gsi = groupset.groups_.front().id_;
  const auto gsf = groupset.groups_.back().id_;

  // Start integration
  std::vector<double> local_leakage(num_gs_groups, 0.0);
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];

    unsigned int f = 0;
    for (const auto& face : cell.faces_)
    {
      if (not face.has_neighbor_ and face.neighbor_id_ == boundary_id)
      {
        const auto& int_f_shape_i = fe_values.intS_shapeI[f];
        const auto num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
        {
          const auto i = cell_mapping.MapFaceNode(f, fi);
          for (unsigned int n = 0; n < num_gs_angles; ++n)
          {
            const auto& omega = quadrature->omegas_[n];
            const auto& weight = quadrature->weights_[n];
            const auto mu = omega.Dot(face.normal_);
            if (mu > 0.0)
            {
              for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
              {
                const auto g = gsg + gsi;
                const auto imap = sdm.MapDOFLocal(cell, i, psi_uk_man, n, g);
                const auto psi = psi_new_local_[groupset_id][imap];
                local_leakage[gsg] += weight * mu * psi * int_f_shape_i[i];
              } // for g
            }   // outgoing
          }     // for n
        }       // for face node
      }         // if right bndry
      ++f;
    } // for face
  }   // for cell

  // Communicate to obtain global leakage
  std::vector<double> global_leakage(num_gs_groups, 0.0);
  mpi_comm.all_reduce(
    local_leakage.data(), num_gs_groups, global_leakage.data(), mpi::op::sum<double>());

  return global_leakage;
}

std::map<uint64_t, std::vector<double>>
DiscreteOrdinatesSolver::ComputeLeakage(const std::vector<uint64_t>& boundary_ids) const
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::ComputeLeakage");

  // Perform checks
  OpenSnLogicalErrorIf(not options_.save_angular_flux,
                       "The option `save_angular_flux` must be set to `true` in order "
                       "to compute outgoing currents.");

  const auto unique_bids = grid_ptr_->GetDomainUniqueBoundaryIDs();
  for (const auto& bid : boundary_ids)
  {
    const auto it = std::find(unique_bids.begin(), unique_bids.end(), bid);
    OpenSnInvalidArgumentIf(it == unique_bids.end(),
                            "Boundary ID " + std::to_string(bid) + "not found on grid.");
  }

  // Initialize local mapping
  std::map<uint64_t, std::vector<double>> local_leakage;
  for (const auto& bid : boundary_ids)
    local_leakage[bid].assign(num_groups_, 0.0);

  // Go through groupsets
  for (unsigned int gs = 0; gs < groupsets_.size(); ++gs)
  {
    const auto& groupset = groupsets_.at(gs);
    const auto& psi_uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature_;

    const auto num_gs_angles = quadrature->omegas_.size();
    const auto num_gs_groups = groupset.groups_.size();
    const auto first_gs_group = groupset.groups_.front().id_;

    const auto& psi_gs = psi_new_local_[gs];

    // Loop over cells for integration
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_mapping = discretization_->GetCellMapping(cell);
      const auto& fe_values = unit_cell_matrices_.at(cell.local_id_);

      unsigned int f = 0;
      for (const auto& face : cell.faces_)
      {
        // If face is on the specified boundary...
        const auto it = std::find(boundary_ids.begin(), boundary_ids.end(), face.neighbor_id_);
        if (not face.has_neighbor_ and it != boundary_ids.end())
        {
          auto& bndry_leakage = local_leakage[face.neighbor_id_];
          const auto& int_f_shape_i = fe_values.intS_shapeI[f];
          const auto num_face_nodes = cell_mapping.NumFaceNodes(f);
          for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
          {
            const auto i = cell_mapping.MapFaceNode(f, fi);
            for (unsigned int n = 0; n < num_gs_angles; ++n)
            {
              const auto& omega = quadrature->omegas_[n];
              const auto& weight = quadrature->weights_[n];
              const auto mu = omega.Dot(face.normal_);
              if (mu <= 0.0)
                continue;

              const auto coeff = weight * mu * int_f_shape_i[i];
              for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
              {
                const auto g = first_gs_group + gsg;
                const auto imap = discretization_->MapDOFLocal(cell, i, psi_uk_man, n, gsg);
                bndry_leakage[g] += coeff * psi_gs[imap];
              } // for groupset group gsg
            }   // for angle n
          }     // for face index fi
        }       // if face on desired boundary
        ++f;
      } // for face
    }   // for cell
  }     // for groupset gs

  // Serialize the data
  std::vector<double> local_data;
  for (const auto& [bid, bndry_leakage] : local_leakage)
    for (const auto& val : bndry_leakage)
      local_data.emplace_back(val);

  // Communicate the data
  std::vector<double> global_data(local_data.size());
  mpi_comm.all_reduce(
    local_data.data(), local_data.size(), global_data.data(), mpi::op::sum<double>());

  // Unpack the data
  std::map<uint64_t, std::vector<double>> global_leakage;
  for (unsigned int b = 0; b < boundary_ids.size(); ++b)
    for (unsigned int g = 0; g < num_groups_; ++g)
      global_leakage[boundary_ids[b]].push_back(global_data[b * num_groups_ + g]);
  return global_leakage;
}

void
DiscreteOrdinatesSolver::InitializeSweepDataStructures()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::InitializeSweepDataStructures");

  log.Log() << program_timer.GetTimeString() << " Initializing sweep datastructures.\n";

  // Define sweep ordering groups
  quadrature_unq_so_grouping_map_.clear();
  std::map<std::shared_ptr<AngularQuadrature>, bool> quadrature_allow_cycles_map_;
  for (auto& groupset : groupsets_)
  {
    if (quadrature_unq_so_grouping_map_.count(groupset.quadrature_) == 0)
      quadrature_unq_so_grouping_map_[groupset.quadrature_] = AssociateSOsAndDirections(
        *grid_ptr_, *groupset.quadrature_, groupset.angleagg_method_, options_.geometry_type);

    if (quadrature_allow_cycles_map_.count(groupset.quadrature_) == 0)
      quadrature_allow_cycles_map_[groupset.quadrature_] = groupset.allow_cycles_;
  }

  // Build sweep orderings
  quadrature_spds_map_.clear();
  for (const auto& [quadrature, info] : quadrature_unq_so_grouping_map_)
  {
    const auto& unique_so_groupings = info.first;

    for (const auto& so_grouping : unique_so_groupings)
    {
      if (so_grouping.empty())
        continue;

      const size_t master_dir_id = so_grouping.front();
      const auto& omega = quadrature->omegas_[master_dir_id];

      bool verbose = false;
      if (not verbose_sweep_angles_.empty())
        for (const size_t dir_id : verbose_sweep_angles_)
          if (VectorListHas(so_grouping, dir_id))
          {
            verbose = true;
            break;
          }

      if (sweep_type_ == "AAH")
      {
        const auto new_swp_order = std::make_shared<SPDS_AdamsAdamsHawkins>(
          omega, *this->grid_ptr_, quadrature_allow_cycles_map_[quadrature], verbose);
        quadrature_spds_map_[quadrature].push_back(new_swp_order);
      }
      else if (sweep_type_ == "CBC")
      {
        const auto new_swp_order = std::make_shared<CBC_SPDS>(
          omega, *this->grid_ptr_, quadrature_allow_cycles_map_[quadrature], verbose);
        quadrature_spds_map_[quadrature].push_back(new_swp_order);
      }
      else
        OpenSnInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
    }
  } // quadrature info-pack

  // Build FLUDS templates
  quadrature_fluds_commondata_map_.clear();
  for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
  {
    for (const auto& spds : spds_list)
    {
      if (sweep_type_ == "AAH")
      {
        quadrature_fluds_commondata_map_[quadrature].push_back(
          std::make_unique<AAH_FLUDSCommonData>(
            grid_nodal_mappings_, *spds, *grid_face_histogram_));
      }
      else if (sweep_type_ == "CBC")
      {
        quadrature_fluds_commondata_map_[quadrature].push_back(
          std::make_unique<CBC_FLUDSCommonData>(*spds, grid_nodal_mappings_));
      }
      else
        OpenSnInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
    }
  } // for quadrature spds-list pair

  log.Log() << program_timer.GetTimeString() << " Done initializing sweep datastructures.\n";
}

std::pair<UniqueSOGroupings, DirIDToSOMap>
DiscreteOrdinatesSolver::AssociateSOsAndDirections(const MeshContinuum& grid,
                                                   const AngularQuadrature& quadrature,
                                                   const AngleAggregationType agg_type,
                                                   const GeometryType lbs_geo_type)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::AssociateSOsAndDirections");

  const std::string fname = __FUNCTION__;

  // Checks
  if (quadrature.omegas_.empty())
    throw std::logic_error(fname + ": Quadrature with no omegas cannot be used.");
  if (quadrature.weights_.empty())
    throw std::logic_error(fname + ": Quadrature with no weights cannot be used.");

  // Build groupings
  UniqueSOGroupings unq_so_grps;
  switch (agg_type)
  {
    // Single
    // The easiest aggregation type. Every direction
    // either has/is assumed to have a unique sweep
    // ordering. Hence there is only group holding ALL
    // the direction indices.
    case AngleAggregationType::SINGLE:
    {
      const size_t num_dirs = quadrature.omegas_.size();
      for (size_t n = 0; n < num_dirs; ++n)
        unq_so_grps.push_back({n});
      break;
    } // case agg_type SINGLE

      // Polar
      // The following conditions allow for polar
      // angle aggregation.
    case AngleAggregationType::POLAR:
    {
      // Check geometry types
      if (not(grid.Type() == ORTHOGONAL or grid.Dimension() == 2 or grid.Extruded()))
        throw std::logic_error(
          fname + ": The simulation is using polar angle aggregation for which only certain "
                  "geometry types are supported, i.e., ORTHOGONAL, 2D or 3D EXTRUDED.");

      // Check quadrature type
      const auto quad_type = quadrature.type_;
      if (quad_type != AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(fname + ": The simulation is using polar angle aggregation for "
                                       "which only Product-type quadratures are supported.");

      // Process Product Quadrature
      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        const auto num_azi = product_quad.azimu_ang_.size();
        const auto num_pol = product_quad.polar_ang_.size();

        // Make two separate list of polar angles
        // One upward-pointing and one downward
        std::vector<size_t> upward_polar_ids;
        std::vector<size_t> dnward_polar_ids;
        for (size_t p = 0; p < num_pol; ++p)
          if (product_quad.polar_ang_[p] > M_PI_2)
            upward_polar_ids.push_back(p);
          else
            dnward_polar_ids.push_back(p);

        // Define lambda working for both upward and dnward polar-ids
        /**Lambda to convert indices and push it onto unq_so_grps.*/
        auto MapPolarAndAzimuthalIDs =
          [&product_quad, &unq_so_grps](const DirIDs& polar_ids, const size_t azimuthal_id)
        {
          DirIDs dir_ids;
          dir_ids.reserve(polar_ids.size());
          for (const size_t p : polar_ids)
            dir_ids.push_back(product_quad.GetAngleNum(p, azimuthal_id));
          unq_so_grps.push_back(std::move(dir_ids));
        };

        // Stack id's for all azimuthal angles
        for (size_t a = 0; a < num_azi; ++a)
        {
          if (not upward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(upward_polar_ids, a);
          if (not dnward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(dnward_polar_ids, a);
        } // for azi-id a

      } // try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(
          fname + ": Casting the angular quadrature to the product quadrature base, failed.");
      }

      break;
    } // case agg_type POLAR

      // Azimuthal
    case AngleAggregationType::AZIMUTHAL:
    {
      // Check geometry types
      if (not(lbs_geo_type == GeometryType::ONED_SPHERICAL or
              lbs_geo_type == GeometryType::TWOD_CYLINDRICAL))
        throw std::logic_error(
          fname + ": The simulation is using azimuthal angle aggregation for which only "
                  "ONED_SPHERICAL or TWOD_CYLINDRICAL derived geometry types are supported.");

      // Check quadrature type
      const auto quad_type = quadrature.type_;
      if (quad_type != AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(fname + ": The simulation is using azimuthal angle aggregation for "
                                       "which only Product-type quadratures are supported.");

      // Process Product Quadrature
      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        for (const auto& dir_set : product_quad.GetDirectionMap())
        {
          std::vector<unsigned int> group1;
          std::vector<unsigned int> group2;
          for (const auto& dir_id : dir_set.second)
            if (quadrature.abscissae_[dir_id].phi > M_PI_2)
              group1.push_back(dir_id);
            else
              group2.push_back(dir_id);

          DirIDs group1_ids(group1.begin(), group1.end());
          DirIDs group2_ids(group2.begin(), group2.end());

          unq_so_grps.push_back(std::move(group1_ids));
          unq_so_grps.push_back(std::move(group2_ids));
        }
      } // try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(
          fname + ": Casting the angular quadrature to the product quadrature base, failed.");
      }

      break;
    }
    default:
      throw std::invalid_argument(fname + ": Called with UNDEFINED angle "
                                          "aggregation type.");
  } // switch angle aggregation type

  // Map directions to sweep orderings
  DirIDToSOMap dir_id_to_so_map;
  {
    size_t so_grouping_id = 0;
    for (const auto& so_grouping : unq_so_grps)
    {
      for (const size_t dir_id : so_grouping)
        dir_id_to_so_map[dir_id] = so_grouping_id;

      ++so_grouping_id;
    } // for so_grouping
  }   // map scope

  return {unq_so_grps, dir_id_to_so_map};
}

void
DiscreteOrdinatesSolver::InitFluxDataStructures(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::InitFluxDataStructures");

  const auto& quadrature_sweep_info = quadrature_unq_so_grouping_map_[groupset.quadrature_];

  const auto& unique_so_groupings = quadrature_sweep_info.first;
  const auto& dir_id_to_so_map = quadrature_sweep_info.second;

  const size_t gs_num_grps = groupset.groups_.size();
  const size_t gs_num_ss = groupset.grp_subset_infos_.size();

  // Passing the sweep boundaries
  //                                            to the angle aggregation
  groupset.angle_agg_ = std::make_shared<AngleAggregation>(
    sweep_boundaries_, num_groups_, gs_num_ss, groupset.quadrature_, grid_ptr_);

  AngleSetGroup angle_set_group;
  size_t angle_set_id = 0;
  for (const auto& so_grouping : unique_so_groupings)
  {
    const size_t master_dir_id = so_grouping.front();
    const size_t so_id = dir_id_to_so_map.at(master_dir_id);

    const auto& sweep_ordering = quadrature_spds_map_[groupset.quadrature_][so_id];
    const auto& fluds_common_data = *quadrature_fluds_commondata_map_[groupset.quadrature_][so_id];

    // Compute direction subsets
    const auto dir_subsets = MakeSubSets(so_grouping.size(), groupset.master_num_ang_subsets_);

    for (size_t gs_ss = 0; gs_ss < gs_num_ss; gs_ss++)
    {
      const size_t gs_ss_size = groupset.grp_subset_infos_[gs_ss].ss_size;
      for (const auto& dir_ss_info : dir_subsets)
      {
        const auto& dir_ss_begin = dir_ss_info.ss_begin;
        const auto& dir_ss_end = dir_ss_info.ss_end;
        const auto& dir_ss_size = dir_ss_info.ss_size;

        std::vector<size_t> angle_indices(dir_ss_size, 0);
        {
          size_t k = 0;
          for (size_t n = dir_ss_begin; n <= dir_ss_end; ++n)
            angle_indices[k++] = so_grouping[n];
        }

        if (sweep_type_ == "AAH")
        {
          std::shared_ptr<FLUDS> fluds = std::make_shared<AAH_FLUDS>(
            gs_ss_size,
            angle_indices.size(),
            dynamic_cast<const AAH_FLUDSCommonData&>(fluds_common_data));

          auto angle_set = std::make_shared<AAH_AngleSet>(angle_set_id++,
                                                          gs_ss_size,
                                                          gs_ss,
                                                          *sweep_ordering,
                                                          fluds,
                                                          angle_indices,
                                                          sweep_boundaries_,
                                                          options_.max_mpi_message_size,
                                                          *grid_local_comm_set_);

          angle_set_group.AngleSets().push_back(angle_set);
        }
        else if (sweep_type_ == "CBC")
        {
          OpenSnLogicalErrorIf(not options_.save_angular_flux,
                               "When using sweep_type \"CBC\" then "
                               "\"save_angular_flux\" must be true.");
          std::shared_ptr<FLUDS> fluds =
            std::make_shared<CBC_FLUDS>(gs_ss_size,
                                        angle_indices.size(),
                                        dynamic_cast<const CBC_FLUDSCommonData&>(fluds_common_data),
                                        psi_new_local_[groupset.id_],
                                        groupset.psi_uk_man_,
                                        *discretization_);

          auto angle_set = std::make_shared<CBC_AngleSet>(angle_set_id++,
                                                          gs_ss_size,
                                                          *sweep_ordering,
                                                          fluds,
                                                          angle_indices,
                                                          sweep_boundaries_,
                                                          gs_ss,
                                                          *grid_local_comm_set_);

          angle_set_group.AngleSets().push_back(angle_set);
        }
        else
          OpenSnInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
      } // for an_ss
    }   // for gs_ss
  }     // for so_grouping

  groupset.angle_agg_->angle_set_groups.push_back(std::move(angle_set_group));

  if (options_.verbose_inner_iterations)
    log.Log() << program_timer.GetTimeString() << " Initialized angle aggregation.";

  opensn::mpi_comm.barrier();
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesSolver::SetSweepChunk(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesSolver::SetSweepChunk");

  if (sweep_type_ == "AAH")
  {
    auto sweep_chunk = std::make_shared<AahSweepChunk>(*grid_ptr_,
                                                       *discretization_,
                                                       unit_cell_matrices_,
                                                       cell_transport_views_,
                                                       densities_local_,
                                                       phi_new_local_,
                                                       psi_new_local_[groupset.id_],
                                                       q_moments_local_,
                                                       groupset,
                                                       matid_to_xs_map_,
                                                       num_moments_,
                                                       max_cell_dof_count_);

    return sweep_chunk;
  }
  else if (sweep_type_ == "CBC")
  {
    auto sweep_chunk = std::make_shared<CbcSweepChunk>(phi_new_local_,
                                                       psi_new_local_[groupset.id_],
                                                       *grid_ptr_,
                                                       *discretization_,
                                                       unit_cell_matrices_,
                                                       cell_transport_views_,
                                                       densities_local_,
                                                       q_moments_local_,
                                                       groupset,
                                                       matid_to_xs_map_,
                                                       num_moments_,
                                                       max_cell_dof_count_);

    return sweep_chunk;
  }
  else
    OpenSnLogicalError("Unsupported sweep_type_ \"" + sweep_type_ + "\"");
}

} // namespace opensn
