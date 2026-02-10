// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/runtime.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "caliper/cali.h"

namespace opensn
{

double
ComputeFissionProduction(LBSProblem& lbs_problem, const std::vector<double>& phi)
{
  CALI_CXX_MARK_SCOPE("ComputeFissionProduction");

  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  const auto& options = lbs_problem.GetOptions();

  const auto first_grp = 0;
  const auto last_grp = lbs_problem.GetNumGroups() - 1;

  // Loop over local cells
  double local_production = 0.0;
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[cell.local_id];

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& F = xs.GetProductionMatrix();
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    if (not xs.IsFissionable())
      continue;

    // Loop over nodes
    const int num_nodes = transport_view.GetNumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const auto uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (unsigned int g = first_grp; g <= last_grp; ++g)
      {
        const auto& prod = F[g];
        for (unsigned int gp = 0; gp <= last_grp; ++gp)
          local_production += prod[gp] * phi[uk_map + gp] * IntV_ShapeI;

        if (options.use_precursors)
          local_production += nu_delayed_sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
      }
    } // for node
  } // for cell

  // Allreduce global production
  double global_production = 0.0;
  mpi_comm.all_reduce(local_production, global_production, mpi::op::sum<double>());

  return global_production;
}

double
ComputeFissionRate(LBSProblem& lbs_problem, const std::vector<double>& phi)
{
  CALI_CXX_MARK_SCOPE("ComputeFissionRate");

  const auto& grid = lbs_problem.GetGrid();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();

  const auto first_grp = 0;
  const auto last_grp = lbs_problem.GetNumGroups() - 1;

  // Loop over local cells
  double local_fission_rate = 0.0;
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices[cell.local_id];

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& sigma_f = xs.GetSigmaFission();

    // skip non-fissionable material
    if (not xs.IsFissionable())
      continue;

    // Loop over nodes
    const int num_nodes = transport_view.GetNumNodes();
    for (auto i = 0; i < num_nodes; ++i)
    {
      const auto uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (unsigned int g = first_grp; g <= last_grp; ++g)
        local_fission_rate += sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
    } // for node
  } // for cell

  // Allreduce global production
  double global_fission_rate = 0.0;
  mpi_comm.all_reduce(local_fission_rate, global_fission_rate, mpi::op::sum<double>());

  return global_fission_rate;
}

void
ComputePrecursors(LBSProblem& lbs_problem)
{
  CALI_CXX_MARK_SCOPE("ComputePrecursors");

  const size_t J = lbs_problem.GetMaxPrecursorsPerMaterial();

  auto& precursor_new_local = lbs_problem.GetPrecursorsNewLocal();
  precursor_new_local.assign(precursor_new_local.size(), 0.0);

  const auto& grid = lbs_problem.GetGrid();
  const auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  const auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto& phi_new_local = lbs_problem.GetPhiNewLocal();

  // Loop over cells
  for (const auto& cell : grid->local_cells)
  {
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = cell_transport_views[cell.local_id];
    const double cell_volume = transport_view.GetVolume();

    // Obtain xs
    const auto& xs = transport_view.GetXS();
    const auto& precursors = xs.GetPrecursors();
    const auto& nu_delayed_sigma_f = xs.GetNuDelayedSigmaF();

    // Loop over precursors
    for (uint64_t j = 0; j < xs.GetNumPrecursors(); ++j)
    {
      size_t dof = cell.local_id * J + j;
      const auto& precursor = precursors[j];
      const double coeff = precursor.fractional_yield / precursor.decay_constant;

      // Loop over nodes
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        const auto uk_map = transport_view.MapDOF(i, 0, 0);
        const double node_V_fraction = fe_values.intV_shapeI(i) / cell_volume;

        // Loop over groups
        for (unsigned int g = 0; g < lbs_problem.GetNumGroups(); ++g)
          precursor_new_local[dof] +=
            coeff * nu_delayed_sigma_f[g] * phi_new_local[uk_map + g] * node_V_fraction;
      } // for node i
    } // for precursor j
  } // for cell
}

void
ComputeBalance(DiscreteOrdinatesProblem& do_problem, double scaling_factor)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::ComputeBalance");

  opensn::mpi_comm.barrier();

  const auto& grid_ = do_problem.GetGrid();
  const auto& discretization_ = do_problem.GetSpatialDiscretization();
  auto& phi_new_local_ = do_problem.GetPhiNewLocal();
  auto& groupsets_ = do_problem.GetGroupsets();
  auto& q_moments_local_ = do_problem.GetQMomentsLocal();
  auto active_set_source_fn = do_problem.GetActiveSetSourceFunction();
  const auto& cell_transport_views_ = do_problem.GetCellTransportViews();
  const auto& unit_cell_matrices_ = do_problem.GetUnitCellMatrices();
  const auto& sweep_boundaries_ = do_problem.GetSweepBoundaries();
  const auto num_groups_ = do_problem.GetNumGroups();
  const auto time_dependent = do_problem.IsTimeDependent();
  const auto dt = time_dependent ? do_problem.GetTimeStep() : 0.0;
  const auto& psi_new_local_ = do_problem.GetPsiNewLocal();
  const auto& psi_old_local_ = do_problem.GetPsiOldLocal();

  // Get material source
  // This is done using the SetSource routine because it allows a lot of flexibility.
  auto mat_src = phi_new_local_;
  mat_src.assign(mat_src.size(), 0.0);
  for (auto& groupset : groupsets_)
  {
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    active_set_source_fn(groupset,
                         q_moments_local_,
                         phi_new_local_,
                         APPLY_FIXED_SOURCES | APPLY_AGS_FISSION_SOURCES |
                           APPLY_WGS_FISSION_SOURCES);
    LBSVecOps::GSScopedCopyPrimarySTLvectors( // NOLINT(readability-suspicious-call-argument)
      do_problem,
      groupset,
      q_moments_local_,
      mat_src);
  }

  // Compute absorption, material-source and in-flow
  double local_out_flow = 0.0;
  double local_in_flow = 0.0;
  double local_absorption = 0.0;
  double local_production = 0.0;
  double local_initial = 0.0;
  double local_final = 0.0;
  for (const auto& cell : grid_->local_cells)
  {
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    const auto& transport_view = cell_transport_views_[cell.local_id];
    const auto& fe_intgrl_values = unit_cell_matrices_[cell.local_id];
    const size_t num_nodes = transport_view.GetNumNodes();
    const auto& IntV_shapeI = fe_intgrl_values.intV_shapeI;
    const auto& IntS_shapeI = fe_intgrl_values.intS_shapeI;

    // Inflow: This is essentially an integration over all faces, all angles, and all groups. For
    // non-reflective boundaries, only the cosines that are negative are added to the inflow
    // integral. For reflective boundaries, it is expected that, upon convergence, inflow = outflow
    // (within numerical tolerances set by the user).
    for (int f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];

      if (not face.has_neighbor) // Boundary face
      {
        const auto& bndry = sweep_boundaries_.at(face.neighbor_id);

        if (bndry->IsReflecting())
        {
          for (unsigned int g = 0; g < num_groups_; ++g)
            local_in_flow += transport_view.GetOutflow(f, g);
        }
        else
        {
          for (const auto& groupset : groupsets_)
          {
            for (int n = 0; n < groupset.quadrature->omegas.size(); ++n)
            {
              const auto& omega = groupset.quadrature->omegas[n];
              const double wt = groupset.quadrature->weights[n];
              const double mu = omega.Dot(face.normal);

              if (mu < 0.0)
              {
                for (int fi = 0; fi < face.vertex_ids.size(); ++fi)
                {
                  const int i = cell_mapping.MapFaceNode(f, fi);
                  const auto& IntFi_shapeI = IntS_shapeI[f](i);

                  for (unsigned int g = groupset.first_group; g <= groupset.last_group; ++g)
                  {
                    const double psi = *bndry->PsiIncoming(cell.local_id, f, fi, n, g);
                    local_in_flow -= mu * wt * psi * IntFi_shapeI;
                  } // for group
                } // for fi
              } // if mu < 0
            } // for n
          } // for groupset
        } // if reflecting boundary
      } // if boundary
    } // for f

    // Outflow: The group-wise outflow was determined during a solve so we just accumulate it here.
    for (size_t f = 0; f < cell.faces.size(); ++f)
      for (unsigned int g = 0; g < num_groups_; ++g)
        local_out_flow += transport_view.GetOutflow(f, g);

    // Absorption and sources
    const auto& xs = transport_view.GetXS();
    const auto& sigma_a = xs.GetSigmaAbsorption();
    const auto& inv_vel = xs.GetInverseVelocity();
    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (unsigned int g = 0; g < num_groups_; ++g)
      {
        auto imap = transport_view.MapDOF(i, 0, g);
        double phi_0g = phi_new_local_[imap];
        double q_0g = mat_src[imap];

        local_absorption += sigma_a[g] * phi_0g * IntV_shapeI(i);
        local_production += q_0g * IntV_shapeI(i);
      } // for g
    } // for i

    if (time_dependent)
    {
      for (const auto& groupset : groupsets_)
      {
        const auto& quad = groupset.quadrature;
        const auto num_angles = quad->omegas.size();
        const auto first_grp = groupset.first_group;
        const auto num_gs_groups = groupset.GetNumGroups();
        const size_t groupset_angle_group_stride =
          groupset.psi_uk_man_.GetNumberOfUnknowns() * num_gs_groups;
        const size_t groupset_group_stride = num_gs_groups;
        const size_t base = discretization_.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0);

        for (size_t i = 0; i < num_nodes; ++i)
        {
          for (size_t gsg = 0; gsg < num_gs_groups; ++gsg)
          {
            double phi_old = 0.0;
            double phi_new = 0.0;
            for (size_t n = 0; n < num_angles; ++n)
            {
              const size_t imap =
                base + i * groupset_angle_group_stride + n * groupset_group_stride;
              const double wt = quad->weights[n];
              phi_old += wt * psi_old_local_[groupset.id][imap + gsg];
              phi_new += wt * psi_new_local_[groupset.id][imap + gsg];
            }

            const auto g = first_grp + gsg;
            const double val = inv_vel[g] * IntV_shapeI(i);
            local_initial += val * phi_old;
            local_final += val * phi_new;
          }
        }
      }
    }
  } // for cell

  // Compute local balance
  double local_balance = local_production + local_in_flow - local_absorption - local_out_flow;
  double local_gain = local_production + local_in_flow;
  std::vector<double> local_balance_table = {
    local_absorption, local_production, local_in_flow, local_out_flow, local_balance, local_gain};
  if (time_dependent)
  {
    local_balance_table.push_back(local_initial);
    local_balance_table.push_back(local_final);
  }
  auto table_size = static_cast<int>(local_balance_table.size());

  // Compute global balance
  std::vector<double> global_balance_table(table_size, 0.0);
  mpi_comm.all_reduce(
    local_balance_table.data(), table_size, global_balance_table.data(), mpi::op::sum<double>());
  double global_absorption = global_balance_table.at(0);
  double global_production = global_balance_table.at(1);
  double global_in_flow = global_balance_table.at(2);
  double global_out_flow = global_balance_table.at(3);
  double global_balance = global_balance_table.at(4);
  double global_gain = global_balance_table.at(5);
  double global_initial = 0.0;
  double global_final = 0.0;
  if (scaling_factor != 1.0)
  {
    global_production *= scaling_factor;
    global_gain = global_production + global_in_flow;
    global_balance = global_gain - (global_absorption + global_out_flow);
  }

  if (time_dependent)
  {
    global_initial = global_balance_table.at(6);
    global_final = global_balance_table.at(7);
  }

  if (time_dependent)
  {
    // NOTE: The time-dependent balance calculation uses end-of-step rates and a simple
    // initial/final calculation. For full consistency with the time discretization, it
    // should use theta-centered rates (sources/absorption), time-averaged boundary
    // inflow/outflow, and a correct handling of pulsed sources/boundaries when
    // start/end times fall inside the timestep.
    const double source = dt * global_production;
    const double absorption = dt * global_absorption;
    const double net_outflow = dt * (global_out_flow - global_in_flow);
    const double net_gain = source - absorption - net_outflow;
    const double conservation_error =
      (global_final == 0.0) ? 0.0 : (global_final - (global_initial + net_gain)) / global_final;

    log.Log() << "\nTimestep balance (dt = " << std::setprecision(6) << std::fixed << dt << "):\n"
              << std::setprecision(6) << std::scientific
              << " Source                      = " << source << "\n"
              << " Absorption                  = " << absorption << "\n"
              << " Net out-flow                = " << net_outflow << "\n"
              << " Initial                     = " << global_initial << "\n"
              << " Final                       = " << global_final << "\n"
              << " Conservation error          = " << conservation_error << "\n\n";
  }
  else
  {
    const double conservation_error = (global_gain == 0.0) ? 0.0 : (global_balance / global_gain);
    log.Log() << "\nBalance table:\n"
              << std::setprecision(6) << std::scientific
              << " Absorption rate             = " << global_absorption << "\n"
              << " Production rate             = " << global_production << "\n"
              << " In-flow rate                = " << global_in_flow << "\n"
              << " Out-flow rate               = " << global_out_flow << "\n"
              << " Gain (In-flow + Production) = " << global_gain << "\n"
              << " Balance (Gain - Loss)       = " << global_balance << "\n"
              << " Conservation error          = " << conservation_error << "\n\n";
  }

  opensn::mpi_comm.barrier();
}

static double
GetInflow(const Cell& cell,
          const CellMapping& cell_mapping,
          const std::vector<UnitCellMatrices>& unit_cell_matrices,
          const CellFace& face,
          const std::shared_ptr<AngularQuadrature> quadrature,
          SweepBoundary& boundary,
          unsigned int f,
          unsigned int g)

{
  const auto& fe_vals = unit_cell_matrices.at(cell.local_id);
  const auto& int_f_shape_i = fe_vals.intS_shapeI[f];
  const unsigned int num_face_nodes = cell_mapping.GetNumFaceNodes(f);
  const size_t num_angles = quadrature->omegas.size();

  double inflow = 0.0;
  for (size_t n = 0; n < num_angles; ++n)
  {
    const auto& omega = quadrature->omegas[n];
    const auto& weight = quadrature->weights[n];
    const double mu = omega.Dot(face.normal);

    if (mu >= 0.0)
      continue;

    for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
    {
      const unsigned int i = cell_mapping.MapFaceNode(f, fi);
      const double* psi_ptr = boundary.PsiIncoming(cell.local_id, f, fi, n, g);
      const double psi_inc = (psi_ptr != nullptr) ? *psi_ptr : 0.0;
      inflow += mu * weight * int_f_shape_i(i) * psi_inc;
    }
  }

  return inflow;
}

std::vector<double>
ComputeLeakage(DiscreteOrdinatesProblem& do_problem,
               const unsigned int groupset_id,
               const uint64_t boundary_id)
{
  CALI_CXX_MARK_SCOPE("ComputeLeakage");

  OpenSnInvalidArgumentIf(groupset_id < 0 or groupset_id >= do_problem.GetGroupsets().size(),
                          "Invalid groupset id.");

  const auto& grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& groupset = do_problem.GetGroupsets().at(groupset_id);
  const auto& cell_transport_views = do_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& sweep_boundaries = do_problem.GetSweepBoundaries();
  const auto& quad = groupset.quadrature;
  const auto num_gs_groups = groupset.GetNumGroups();
  const auto gsi = groupset.first_group;

  std::vector<double> local_leakage(num_gs_groups, 0.0);
  for (const auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);

    for (unsigned int f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];

      if (not face.has_neighbor && face.neighbor_id == boundary_id)
      {
        auto& bndry = *sweep_boundaries.at(face.neighbor_id);

        for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
        {
          auto g = gsi + gsg;
          local_leakage[gsg] += transport_view.GetOutflow(f, g);
          local_leakage[gsg] +=
            GetInflow(cell, cell_mapping, unit_cell_matrices, face, quad, bndry, f, g);
        }
      }
    }
  }

  // Reduce
  std::vector<double> global_leakage(num_gs_groups, 0.0);
  mpi_comm.all_reduce(local_leakage.data(),
                      static_cast<int>(num_gs_groups),
                      global_leakage.data(),
                      mpi::op::sum<double>());

  return global_leakage;
}

std::map<uint64_t, std::vector<double>>
ComputeLeakage(DiscreteOrdinatesProblem& do_problem, const std::vector<uint64_t>& boundary_ids)
{
  CALI_CXX_MARK_SCOPE("ComputeLeakage");

  const auto& grid = do_problem.GetGrid();
  const auto uniq = grid->GetUniqueBoundaryIDs();
  for (const auto& bid : boundary_ids)
    OpenSnInvalidArgumentIf(std::find(uniq.begin(), uniq.end(), bid) == uniq.end(),
                            "Boundary id " + std::to_string(bid) + " not found on grid.");

  const auto num_groups = do_problem.GetNumGroups();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();
  const auto& cell_transport_views = do_problem.GetCellTransportViews();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& sweep_boundaries = do_problem.GetSweepBoundaries();

  std::map<uint64_t, std::vector<double>> local_leakage;
  for (const auto& bid : boundary_ids)
    local_leakage[bid].assign(num_groups, 0.0);

  const std::unordered_set<uint64_t> requested_ids(boundary_ids.begin(), boundary_ids.end());

  for (const auto& groupset : groupsets)
  {
    const auto& quad = groupset.quadrature;

    for (const auto& cell : grid->local_cells)
    {
      const auto& transport_view = cell_transport_views[cell.local_id];
      const auto& cell_mapping = sdm.GetCellMapping(cell);

      for (unsigned int f = 0; f < cell.faces.size(); ++f)
      {
        const auto& face = cell.faces[f];

        if (not face.has_neighbor && requested_ids.count(face.neighbor_id))
        {
          auto& bndry = *sweep_boundaries.at(face.neighbor_id);

          for (auto g = groupset.first_group; g <= groupset.last_group; ++g)
          {
            local_leakage[face.neighbor_id][g] += transport_view.GetOutflow(f, g);
            local_leakage[face.neighbor_id][g] +=
              GetInflow(cell, cell_mapping, unit_cell_matrices, face, quad, bndry, f, g);
          }
        }
      }
    }
  }

  // Serialize in user-requested order
  std::vector<double> local_data;
  local_data.reserve(boundary_ids.size() * num_groups);
  for (const auto& bid : boundary_ids)
  {
    const auto& vec = local_leakage.at(bid);
    local_data.insert(local_data.end(), vec.begin(), vec.end());
  }

  // Reduce
  auto count = static_cast<const int>(local_data.size());
  std::vector<double> global_data(count, 0.0);
  mpi_comm.all_reduce(local_data.data(), count, global_data.data(), mpi::op::sum<double>());

  // Unpack in user-requested order
  std::map<uint64_t, std::vector<double>> global_leakage;
  for (const auto& bid : boundary_ids)
    global_leakage[bid] = std::vector<double>();

  for (size_t b = 0; b < boundary_ids.size(); ++b)
  {
    auto& vec = global_leakage[boundary_ids[b]];
    vec.reserve(num_groups);
    for (unsigned int g = 0; g < num_groups; ++g)
      vec.push_back(global_data[b * num_groups + g]);
  }

  return global_leakage;
}

} // namespace opensn
