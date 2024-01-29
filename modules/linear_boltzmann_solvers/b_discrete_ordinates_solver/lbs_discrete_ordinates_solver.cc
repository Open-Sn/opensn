#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "framework/object_factory.h"
#include "framework/memory_usage.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "framework/mesh/sweep_utilities/spds/spds_adams_adams_hawkins.h"
#include "framework/mesh/sweep_utilities/fluds/aah_fluds.h"
#include "framework/mesh/sweep_utilities/angle_set/aah_angle_set.h"
#include "framework/math/quadratures/angular_quadrature_base.h"
#include "framework/math/quadratures/angular_product_quadrature.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/source_functions/source_function.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_structs.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_spds.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep_chunks/cbc_sweep_chunk.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/logging/log_exceptions.h"
#include <iomanip>

namespace opensn
{
namespace lbs
{

OpenSnRegisterObject(lbs, DiscreteOrdinatesSolver);

DiscreteOrdinatesSolver::DiscreteOrdinatesSolver(const std::string& text_name)
  : LBSSolver(text_name)
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
  for (auto& groupset : groupsets_)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);
  }
}

std::pair<size_t, size_t>
DiscreteOrdinatesSolver::GetNumPhiIterativeUnknowns()
{
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

  source_event_tag_ = log.GetRepeatingEventTag("Set Source");
}

void
DiscreteOrdinatesSolver::InitializeWGSSolvers()
{
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

    auto wgs_solver = std::make_shared<WGSLinearSolver>(sweep_wgs_context_ptr);

    wgs_solvers_.push_back(wgs_solver);
  } // for groupset
}

void
DiscreteOrdinatesSolver::ScalePhiVector(PhiSTLOption which_phi, double value)
{
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
DiscreteOrdinatesSolver::SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset,
                                                           Vec x,
                                                           PhiSTLOption which_phi)
{
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
DiscreteOrdinatesSolver::SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset,
                                                           Vec x_src,
                                                           PhiSTLOption which_phi)
{
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
  VecGetArrayRead(x_src, &x_ref);

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

  VecRestoreArrayRead(x_src, &x_ref);
}

void
DiscreteOrdinatesSolver::GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                                       PhiSTLOption from_which_phi,
                                                       PhiSTLOption to_which_phi)
{
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
DiscreteOrdinatesSolver::SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& gs_ids,
                                                                Vec x,
                                                                PhiSTLOption which_phi)
{
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
  for (int gs_id : gs_ids)
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
DiscreteOrdinatesSolver::SetPrimarySTLvectorFromMultiGSPETScVecFrom(const std::vector<int>& gs_ids,
                                                                    Vec x_src,
                                                                    PhiSTLOption which_phi)
{
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
  VecGetArrayRead(x_src, &x_ref);

  int64_t index = -1;
  for (int gs_id : gs_ids)
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

  VecRestoreArrayRead(x_src, &x_ref);
}

void
DiscreteOrdinatesSolver::ZeroOutflowBalanceVars(LBSGroupset& groupset)
{
  for (auto& cell_transport_view : cell_transport_views_)
    for (auto& group : groupset.groups_)
      cell_transport_view.ZeroOutflow(group.id_);
}

void
DiscreteOrdinatesSolver::ComputeBalance()
{
  opensn::mpi.Barrier();
  log.Log() << "\n********** Computing balance\n";

  // Get material source
  // This is done using the SetSource routine
  // because it allows a lot of flexibility.
  auto mat_src = phi_old_local_;
  mat_src.assign(mat_src.size(), 0.0);
  for (auto& groupset : groupsets_)
  {
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    active_set_source_function_(groupset,
                                q_moments_local_,
                                PhiOldLocal(),
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

    // Inflow
    // This is essentially an integration over
    // all faces, all angles, and all groups.
    // Only the cosines that are negative are
    // added to the integral.
    for (int f = 0; f < cell.faces_.size(); ++f)
    {
      const auto& face = cell.faces_[f];

      for (const auto& groupset : groupsets_)
      {
        for (int n = 0; n < groupset.quadrature_->omegas_.size(); ++n)
        {
          const auto& omega = groupset.quadrature_->omegas_[n];
          const double wt = groupset.quadrature_->weights_[n];
          const double mu = omega.Dot(face.normal_);

          if (mu < 0.0 and (not face.has_neighbor_)) // mu<0 and bndry
          {
            const auto& bndry = sweep_boundaries_[face.neighbor_id_];
            for (int fi = 0; fi < face.vertex_ids_.size(); ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);
              const auto& IntFi_shapeI = IntS_shapeI[f][i];

              for (const auto& group : groupset.groups_)
              {
                const int g = group.id_;
                const double psi = *bndry->HeterogeneousPsiIncoming(cell.local_id_, f, fi, n, g, 0);
                local_in_flow -= mu * wt * psi * IntFi_shapeI;
              } // for g
            }   // for fi
          }     // if bndry
        }       // for n
      }         // for groupset
    }           // for f

    // Outflow
    // The group-wise outflow was determined
    // during a solve so here we just
    // consolidate it.
    for (int g = 0; g < num_groups_; ++g)
      local_out_flow += transport_view.GetOutflow(g);

    // Absorption and Src
    // Isotropic flux based absorption and source
    const auto& xs = transport_view.XS();
    const auto& sigma_a = xs.SigmaAbsorption();
    for (int i = 0; i < num_nodes; ++i)
      for (int g = 0; g < num_groups_; ++g)
      {
        size_t imap = transport_view.MapDOF(i, 0, g);
        double phi_0g = phi_old_local_[imap];
        double q_0g = mat_src[imap];

        local_absorption += sigma_a[g] * phi_0g * IntV_shapeI[i];
        local_production += q_0g * IntV_shapeI[i];
      } // for g
  }     // for cell

  // Consolidate local balances
  double local_balance = local_production + local_in_flow - local_absorption - local_out_flow;
  double local_gain = local_production + local_in_flow;

  std::vector<double> local_balance_table = {
    local_absorption, local_production, local_in_flow, local_out_flow, local_balance, local_gain};
  size_t table_size = local_balance_table.size();

  std::vector<double> globl_balance_table(table_size, 0.0);

  MPI_Allreduce(local_balance_table.data(),
                globl_balance_table.data(),
                table_size,
                MPI_DOUBLE,
                MPI_SUM,
                mpi.comm);

  double globl_absorption = globl_balance_table.at(0);
  double globl_production = globl_balance_table.at(1);
  double globl_in_flow = globl_balance_table.at(2);
  double globl_out_flow = globl_balance_table.at(3);
  double globl_balance = globl_balance_table.at(4);
  double globl_gain = globl_balance_table.at(5);

  log.Log() << "Balance table:\n"
            << std::setprecision(5) << std::scientific
            << " Absorption rate              = " << globl_absorption << "\n"
            << " Production rate              = " << globl_production << "\n"
            << " In-flow rate                 = " << globl_in_flow << "\n"
            << " Out-flow rate                = " << globl_out_flow << "\n"
            << " Net Gain (In-flow + sources) = " << globl_gain << "\n"
            << " Net Balance                  = " << globl_balance << "\n"
            << " (Net Balance)/(Net Gain)     = " << globl_balance / globl_gain << "\n";

  log.Log() << "\n********** Done computing balance\n";

  opensn::mpi.Barrier();
}

std::vector<double>
DiscreteOrdinatesSolver::ComputeLeakage(const unsigned int groupset_id,
                                        const uint64_t boundary_id) const
{
  // Perform checks
  ChiInvalidArgumentIf(groupset_id < 0 or groupset_id >= groupsets_.size(), "Invalid groupset id.");
  ChiLogicalErrorIf(not options_.save_angular_flux,
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
  // clang-format off
  std::vector<double> global_leakage(num_gs_groups, 0.0);
  MPI_Allreduce(local_leakage.data(),   // sendbuf
                global_leakage.data(),  // recvbuf
                num_gs_groups,          // count
                MPI_DOUBLE,             // type
                MPI_SUM,                // operation
                opensn::mpi.comm);         // comm
  // clang-format on

  return global_leakage;
}

std::map<uint64_t, std::vector<double>>
DiscreteOrdinatesSolver::ComputeLeakage(const std::vector<uint64_t>& boundary_ids) const
{
  // Perform checks
  ChiLogicalErrorIf(not options_.save_angular_flux,
                    "The option `save_angular_flux` must be set to `true` in order "
                    "to compute outgoing currents.");

  const auto unique_bids = grid_ptr_->GetDomainUniqueBoundaryIDs();
  for (const auto& bid : boundary_ids)
  {
    const auto it = std::find(unique_bids.begin(), unique_bids.end(), bid);
    ChiInvalidArgumentIf(it == unique_bids.end(),
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
              if (mu <= 0.0) continue;

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
  // clang-format off
  std::vector<double> global_data(local_data.size());
  MPI_Allreduce(local_data.data(),    // sendbuf
                global_data.data(),   // recvbuf
                local_data.size(),    // count
                MPI_DOUBLE,           // type
                MPI_SUM,              // operation
                opensn::mpi.comm);       // comm
  // clang-format on

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
  log.Log() << program_timer.GetTimeString() << " Initializing sweep datastructures.\n";

  // Perform checks
  {
    auto& mesh_handler = GetCurrentHandler();
    auto& mesher = mesh_handler.GetVolumeMesher();

    for (const auto& groupset : groupsets_)
    {
      bool no_cycles_parmetis_partitioning =
        ((mesher.options.partition_type == VolumeMesher::PartitionType::PARMETIS) and
         (not groupset.allow_cycles_));

      bool is_1D_geometry = options_.geometry_type == GeometryType::ONED_SLAB;

      if (no_cycles_parmetis_partitioning and not is_1D_geometry and
          (opensn::mpi.process_count > 1))
        throw std::logic_error("When using PARMETIS type partitioning then groupset iterative "
                               "method must be NPT_CLASSICRICHARDSON_CYCLES or NPT_GMRES_CYCLES");
    } // for groupset
  }

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
      if (so_grouping.empty()) continue;

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
        ChiInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
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
        ChiInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
    }
  } // for quadrature spds-list pair

  log.Log() << program_timer.GetTimeString() << " Done initializing sweep datastructures.\n";
}

std::pair<UniqueSOGroupings, DirIDToSOMap>
DiscreteOrdinatesSolver::AssociateSOsAndDirections(const MeshContinuum& grid,
                                                   const AngularQuadrature& quadrature,
                                                   const AngleAggregationType agg_type,
                                                   const lbs::GeometryType lbs_geo_type)
{
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
      const auto grid_attribs = grid.Attributes();
      if (not(grid_attribs & ORTHOGONAL or grid_attribs & DIMENSION_2 or grid_attribs & EXTRUDED))
        throw std::logic_error(
          fname + ": The simulation is using polar angle aggregation for which only certain "
                  "geometry types are supported, i.e., ORTHOGONAL, DIMENSION_2 or 3D EXTRUDED.");

      // Check quadrature type
      const auto quad_type = quadrature.type_;
      if (quad_type != AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(fname + ": The simulation is using polar angle aggregation for "
                                       "which only Product-type quadratures are supported.");

      // Process Product Quadrature
      try
      {
        typedef ProductQuadrature ProdQuadType;
        const auto& product_quad = dynamic_cast<const ProdQuadType&>(quadrature);

        const auto num_azi = product_quad.azimu_ang_.size();
        const auto num_pol = product_quad.polar_ang_.size();

        // Make two separate list of polar angles
        // One upward-pointing and one downward
        std::vector<size_t> upward_polar_ids;
        std::vector<size_t> dnward_polar_ids;
        for (size_t p = 0; p < num_pol; ++p)
          if (product_quad.polar_ang_[p] > M_PI_2) upward_polar_ids.push_back(p);
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
          if (not upward_polar_ids.empty()) MapPolarAndAzimuthalIDs(upward_polar_ids, a);
          if (not dnward_polar_ids.empty()) MapPolarAndAzimuthalIDs(dnward_polar_ids, a);
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
        typedef ProductQuadrature ProdQuadType;
        const auto& product_quad = dynamic_cast<const ProdQuadType&>(quadrature);

        for (const auto& dir_set : product_quad.GetDirectionMap())
        {
          std::vector<unsigned int> group1;
          std::vector<unsigned int> group2;
          for (const auto& dir_id : dir_set.second)
            if (quadrature.abscissae_[dir_id].phi > M_PI_2) group1.push_back(dir_id);
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
  typedef AngleSetGroup TAngleSetGroup;
  typedef AAH_AngleSet TAAH_AngleSet;

  const auto& quadrature_sweep_info = quadrature_unq_so_grouping_map_[groupset.quadrature_];

  const auto& unique_so_groupings = quadrature_sweep_info.first;
  const auto& dir_id_to_so_map = quadrature_sweep_info.second;

  const size_t gs_num_grps = groupset.groups_.size();
  const size_t gs_num_ss = groupset.grp_subset_infos_.size();

  // Passing the sweep boundaries
  //                                            to the angle aggregation
  typedef AngleAggregation AngleAgg;
  groupset.angle_agg_ = std::make_shared<AngleAgg>(
    sweep_boundaries_, gs_num_grps, gs_num_ss, groupset.quadrature_, grid_ptr_);

  TAngleSetGroup angle_set_group;
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

          auto angleSet = std::make_shared<TAAH_AngleSet>(angle_set_id++,
                                                          gs_ss_size,
                                                          gs_ss,
                                                          *sweep_ordering,
                                                          fluds,
                                                          angle_indices,
                                                          sweep_boundaries_,
                                                          options_.sweep_eager_limit,
                                                          *grid_local_comm_set_);

          angle_set_group.AngleSets().push_back(angleSet);
        }
        else if (sweep_type_ == "CBC")
        {
          ChiLogicalErrorIf(not options_.save_angular_flux,
                            "When using sweep_type \"CBC\" then "
                            "\"save_angular_flux\" must be true.");
          std::shared_ptr<FLUDS> fluds =
            std::make_shared<CBC_FLUDS>(gs_ss_size,
                                        angle_indices.size(),
                                        dynamic_cast<const CBC_FLUDSCommonData&>(fluds_common_data),
                                        psi_new_local_[groupset.id_],
                                        groupset.psi_uk_man_,
                                        *discretization_);

          auto angleSet = std::make_shared<CBC_AngleSet>(angle_set_id++,
                                                         gs_ss_size,
                                                         *sweep_ordering,
                                                         fluds,
                                                         angle_indices,
                                                         sweep_boundaries_,
                                                         gs_ss,
                                                         *grid_local_comm_set_);

          angle_set_group.AngleSets().push_back(angleSet);
        }
        else
          ChiInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
      } // for an_ss
    }   // for gs_ss
  }     // for so_grouping

  groupset.angle_agg_->angle_set_groups.push_back(std::move(angle_set_group));

  if (options_.verbose_inner_iterations)
    log.Log() << program_timer.GetTimeString() << " Initialized Angle Aggregation.   "
              << "         Process memory = " << std::setprecision(3) << GetMemoryUsageInMB()
              << " MB.";

  opensn::mpi.Barrier();
}

void
DiscreteOrdinatesSolver::ResetSweepOrderings(LBSGroupset& groupset)
{
  log.Log0Verbose1() << "Resetting SPDS and FLUDS";

  groupset.angle_agg_->angle_set_groups.clear();

  opensn::mpi.Barrier();

  log.Log() << "SPDS and FLUDS reset complete.            Process memory = " << std::setprecision(3)
            << GetMemoryUsageInMB() << " MB";

  double local_app_memory =
    log.ProcessEvent(Logger::StdTags::MAX_MEMORY_USAGE, Logger::EventOperation::MAX_VALUE);
  double total_app_memory = 0.0;
  MPI_Allreduce(&local_app_memory, &total_app_memory, 1, MPI_DOUBLE, MPI_SUM, mpi.comm);
  double max_proc_memory = 0.0;
  MPI_Allreduce(&local_app_memory, &max_proc_memory, 1, MPI_DOUBLE, MPI_MAX, mpi.comm);

  log.Log() << "\n"
            << std::setprecision(3)
            << "           Total application memory (max): " << total_app_memory / 1024.0 << " GB\n"
            << "           Maximum process memory        : " << max_proc_memory / 1024.0
            << " GB\n\n";
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesSolver::SetSweepChunk(LBSGroupset& groupset)
{
  if (sweep_type_ == "AAH")
  {
    auto sweep_chunk = std::make_shared<AAH_SweepChunk>(*grid_ptr_,
                                                        *discretization_,
                                                        unit_cell_matrices_,
                                                        cell_transport_views_,
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
    auto sweep_chunk = std::make_shared<CBC_SweepChunk>(phi_new_local_,
                                                        psi_new_local_[groupset.id_],
                                                        *grid_ptr_,
                                                        *discretization_,
                                                        unit_cell_matrices_,
                                                        cell_transport_views_,
                                                        q_moments_local_,
                                                        groupset,
                                                        matid_to_xs_map_,
                                                        num_moments_,
                                                        max_cell_dof_count_);

    return sweep_chunk;
  }
  else
    ChiLogicalError("Unsupported sweep_type_ \"" + sweep_type_ + "\"");
}

} // namespace lbs
} // namespace opensn
