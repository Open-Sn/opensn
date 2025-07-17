// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "caliper/cali.h"

namespace opensn
{

double
ComputeFissionProduction(LBSProblem& lbs_problem, const std::vector<double>& phi)
{
  CALI_CXX_MARK_SCOPE("ComputeFissionProduction");

  auto& groups = lbs_problem.GetGroups();
  auto& grid = lbs_problem.GetGrid();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  auto& options = lbs_problem.GetOptions();

  const int first_grp = groups.front().id;
  const int last_grp = groups.back().id;

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
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
      {
        const auto& prod = F[g];
        for (size_t gp = 0; gp <= last_grp; ++gp)
          local_production += prod[gp] * phi[uk_map + gp] * IntV_ShapeI;

        if (options.use_precursors)
          for (unsigned int j = 0; j < xs.GetNumPrecursors(); ++j)
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

  auto& groups = lbs_problem.GetGroups();
  auto& grid = lbs_problem.GetGrid();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
  auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();

  const int first_grp = groups.front().id;
  const int last_grp = groups.back().id;

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
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.intV_shapeI(i);

      // Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
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

  auto& grid = lbs_problem.GetGrid();
  auto& groups = lbs_problem.GetGroups();
  auto& unit_cell_matrices = lbs_problem.GetUnitCellMatrices();
  auto& cell_transport_views = lbs_problem.GetCellTransportViews();
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
        const size_t uk_map = transport_view.MapDOF(i, 0, 0);
        const double node_V_fraction = fe_values.intV_shapeI(i) / cell_volume;

        // Loop over groups
        for (unsigned int g = 0; g < groups.size(); ++g)
          precursor_new_local[dof] +=
            coeff * nu_delayed_sigma_f[g] * phi_new_local[uk_map + g] * node_V_fraction;
      } // for node i
    } // for precursor j
  } // for cell
}

void
ComputeBalance(DiscreteOrdinatesProblem& do_problem)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::ComputeBalance");

  opensn::mpi_comm.barrier();

  auto& grid_ = do_problem.GetGrid();
  auto& discretization_ = do_problem.GetSpatialDiscretization();
  auto& phi_new_local_ = do_problem.GetPhiNewLocal();
  auto& groupsets_ = do_problem.GetGroupsets();
  auto& q_moments_local_ = do_problem.GetQMomentsLocal();
  auto active_set_source_fn = do_problem.GetActiveSetSourceFunction();
  auto& cell_transport_views_ = do_problem.GetCellTransportViews();
  auto& unit_cell_matrices_ = do_problem.GetUnitCellMatrices();
  auto& sweep_boundaries_ = do_problem.GetSweepBoundaries();
  const auto num_groups_ = do_problem.GetNumGroups();

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
    LBSVecOps::GSScopedCopyPrimarySTLvectors(do_problem, groupset, q_moments_local_, mat_src);
  }

  // Compute absorption, material-source and in-flow
  double local_out_flow = 0.0;
  double local_in_flow = 0.0;
  double local_absorption = 0.0;
  double local_production = 0.0;
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
          for (int g = 0; g < num_groups_; ++g)
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

                  for (const auto& group : groupset.groups)
                  {
                    const int g = group.id;
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
    for (int f = 0; f < cell.faces.size(); ++f)
      for (int g = 0; g < num_groups_; ++g)
        local_out_flow += transport_view.GetOutflow(f, g);

    // Absorption and sources
    const auto& xs = transport_view.GetXS();
    const auto& sigma_a = xs.GetSigmaAbsorption();
    for (int i = 0; i < num_nodes; ++i)
    {
      for (int g = 0; g < num_groups_; ++g)
      {
        size_t imap = transport_view.MapDOF(i, 0, g);
        double phi_0g = phi_new_local_[imap];
        double q_0g = mat_src[imap];

        local_absorption += sigma_a[g] * phi_0g * IntV_shapeI(i);
        local_production += q_0g * IntV_shapeI(i);
      } // for g
    } // for i
  } // for cell

  // Compute local balance
  double local_balance = local_production + local_in_flow - local_absorption - local_out_flow;
  double local_gain = local_production + local_in_flow;
  std::vector<double> local_balance_table = {
    local_absorption, local_production, local_in_flow, local_out_flow, local_balance, local_gain};
  size_t table_size = local_balance_table.size();

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

  log.Log() << "Balance table:\n"
            << std::setprecision(6) << std::scientific
            << " Absorption rate             = " << global_absorption << "\n"
            << " Production rate             = " << global_production << "\n"
            << " In-flow rate                = " << global_in_flow << "\n"
            << " Out-flow rate               = " << global_out_flow << "\n"
            << " Gain (In-flow + Production) = " << global_gain << "\n"
            << " Balance (Gain - Loss)       = " << global_balance << "\n"
            << " Balance/Gain, in %          = " << global_balance / global_gain * 100. << "\n";

  opensn::mpi_comm.barrier();
}

std::vector<double>
ComputeLeakage(DiscreteOrdinatesProblem& do_problem,
               const unsigned int groupset_id,
               const uint64_t boundary_id)
{
  CALI_CXX_MARK_SCOPE("ComputeLeakage");

  // Perform checks
  OpenSnInvalidArgumentIf(groupset_id < 0 or groupset_id >= do_problem.GetGroupsets().size(),
                          "Invalid groupset id.");
  OpenSnLogicalErrorIf(not do_problem.GetOptions().save_angular_flux,
                       "The option `save_angular_flux` must be set to `true` in order "
                       "to compute outgoing currents.");

  const auto& grid = do_problem.GetGrid();
  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& groupset = do_problem.GetGroupsets().at(groupset_id);
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& psi_new_local = do_problem.GetPsiNewLocal();
  const auto& psi_uk_man = groupset.psi_uk_man_;
  const auto& quadrature = groupset.quadrature;

  const auto num_gs_angles = quadrature->omegas.size();
  const auto num_gs_groups = groupset.groups.size();

  const auto gsi = groupset.groups.front().id;
  const auto gsf = groupset.groups.back().id;

  // Start integration
  std::vector<double> local_leakage(num_gs_groups, 0.0);
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& fe_values = unit_cell_matrices[cell.local_id];

    unsigned int f = 0;
    for (const auto& face : cell.faces)
    {
      if (not face.has_neighbor and face.neighbor_id == boundary_id)
      {
        const auto& int_f_shape_i = fe_values.intS_shapeI[f];
        const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
        {
          const auto i = cell_mapping.MapFaceNode(f, fi);
          for (unsigned int n = 0; n < num_gs_angles; ++n)
          {
            const auto& omega = quadrature->omegas[n];
            const auto& weight = quadrature->weights[n];
            const auto mu = omega.Dot(face.normal);
            if (mu > 0.0)
            {
              for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
              {
                const auto g = gsg + gsi;
                const auto imap = sdm.MapDOFLocal(cell, i, psi_uk_man, n, g);
                const auto psi = psi_new_local[groupset_id][imap];
                local_leakage[gsg] += weight * mu * psi * int_f_shape_i(i);
              } // for g
            } // outgoing
          } // for n
        } // for face node
      } // if right bndry
      ++f;
    } // for face
  } // for cell

  // Communicate to obtain global leakage
  std::vector<double> global_leakage(num_gs_groups, 0.0);
  mpi_comm.all_reduce(
    local_leakage.data(), num_gs_groups, global_leakage.data(), mpi::op::sum<double>());

  return global_leakage;
}

std::map<uint64_t, std::vector<double>>
ComputeLeakage(DiscreteOrdinatesProblem& do_problem, const std::vector<uint64_t>& boundary_ids)
{
  CALI_CXX_MARK_SCOPE("ComputeLeakage");

  // Perform checks
  OpenSnLogicalErrorIf(not do_problem.GetOptions().save_angular_flux,
                       "The option `save_angular_flux` must be set to `true` in order "
                       "to compute outgoing currents.");

  const auto& grid = do_problem.GetGrid();
  const auto unique_bids = grid->GetUniqueBoundaryIDs();
  for (const auto& bid : boundary_ids)
  {
    const auto it = std::find(unique_bids.begin(), unique_bids.end(), bid);
    OpenSnInvalidArgumentIf(it == unique_bids.end(),
                            "Boundary ID " + std::to_string(bid) + "not found on grid.");
  }

  const auto num_groups = do_problem.GetNumGroups();
  // Initialize local mapping
  std::map<uint64_t, std::vector<double>> local_leakage;
  for (const auto& bid : boundary_ids)
    local_leakage[bid].assign(num_groups, 0.0);

  const auto& sdm = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const auto& psi_new_local = do_problem.GetPsiNewLocal();
  // Go through groupsets
  for (unsigned int gs = 0; gs < groupsets.size(); ++gs)
  {
    const auto& groupset = groupsets.at(gs);
    const auto& psi_uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;

    const auto num_gs_angles = quadrature->omegas.size();
    const auto num_gs_groups = groupset.groups.size();
    const auto first_gs_group = groupset.groups.front().id;

    const auto& psi_gs = psi_new_local[gs];

    // Loop over cells for integration
    for (const auto& cell : grid->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const auto& fe_values = unit_cell_matrices.at(cell.local_id);

      unsigned int f = 0;
      for (const auto& face : cell.faces)
      {
        // If face is on the specified boundary...
        const auto it = std::find(boundary_ids.begin(), boundary_ids.end(), face.neighbor_id);
        if (not face.has_neighbor and it != boundary_ids.end())
        {
          auto& bndry_leakage = local_leakage[face.neighbor_id];
          const auto& int_f_shape_i = fe_values.intS_shapeI[f];
          const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
          for (unsigned int fi = 0; fi < num_face_nodes; ++fi)
          {
            const auto i = cell_mapping.MapFaceNode(f, fi);
            for (unsigned int n = 0; n < num_gs_angles; ++n)
            {
              const auto& omega = quadrature->omegas[n];
              const auto& weight = quadrature->weights[n];
              const auto mu = omega.Dot(face.normal);
              if (mu <= 0.0)
                continue;

              const auto coeff = weight * mu * int_f_shape_i(i);
              for (unsigned int gsg = 0; gsg < num_gs_groups; ++gsg)
              {
                const auto g = first_gs_group + gsg;
                const auto imap = sdm.MapDOFLocal(cell, i, psi_uk_man, n, gsg);
                bndry_leakage[g] += coeff * psi_gs[imap];
              } // for groupset group gsg
            } // for angle n
          } // for face index fi
        } // if face on desired boundary
        ++f;
      } // for face
    } // for cell
  } // for groupset gs

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
    for (unsigned int g = 0; g < num_groups; ++g)
      global_leakage[boundary_ids[b]].push_back(global_data[b * num_groups + g]);
  return global_leakage;
}

} // namespace opensn
