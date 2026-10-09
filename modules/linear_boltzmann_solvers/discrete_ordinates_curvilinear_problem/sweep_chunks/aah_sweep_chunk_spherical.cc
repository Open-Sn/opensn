// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/discrete_ordinates_curvilinear_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_curvilinear_problem/sweep_chunks/aah_sweep_chunk_spherical.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <cmath>
#include <stdexcept>

namespace opensn
{

AAHSweepChunkSpherical::AAHSweepChunkSpherical(DiscreteOrdinatesProblem& problem,
                                               LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetCellOutflowViews(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    secondary_unit_cell_matrices_(dynamic_cast<const DiscreteOrdinatesCurvilinearProblem&>(problem)
                                    .GetSecondaryUnitCellMatrices())
{
  if (std::dynamic_pointer_cast<GLProductQuadrature1DSpherical>(groupset_.quadrature) == nullptr)
    throw std::invalid_argument(
      "AAHSweepChunkSpherical: 1D spherical geometry requires a GLProductQuadrature1DSpherical.");

  unknown_manager_.AddUnknown(UnknownType::VECTOR_N, groupset_.GetNumGroups());
  psi_edge_.assign(discretization_.GetNumLocalDOFs(unknown_manager_), 0.0);
}

void
AAHSweepChunkSpherical::Sweep(AngleSet& angle_set)
{
  const auto gs_size = groupset_.GetNumGroups();
  const auto gs_gi = groupset_.first_group;

  int deploc_face_counter = -1;
  int preloc_face_counter = -1;

  auto& fluds = dynamic_cast<AAH_FLUDS&>(angle_set.GetFLUDS());
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();
  const auto quadrature =
    std::dynamic_pointer_cast<CurvilinearProductQuadrature>(groupset_.quadrature);
  const auto& fac_diamond_difference = quadrature->GetDiamondDifferenceFactor();
  const auto& fac_streaming_operator = quadrature->GetStreamingOperatorFactor();

  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(gs_size, Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& spds = angle_set.GetSPDS();
  const auto& spls = spds.GetLocalSubgrid();
  const auto& as_angle_indices = angle_set.GetAngleIndices();
  for (size_t spls_index = 0; spls_index < spls.size(); ++spls_index)
  {
    const auto cell_local_id = spls[spls_index];
    const auto& cell = grid_->GetLocalCell(cell_local_id);
    const auto& cell_mapping = discretization_.GetCellMapping(cell);
    const auto& cell_transport_view = cell_transport_views_[cell_local_id];
    auto& cell_outflow_view = cell_outflow_views_[cell_local_id];
    const auto cell_num_faces = cell.faces.size();
    const auto cell_num_nodes = cell_mapping.GetNumNodes();
    const auto& face_orientations = spds.GetCellFaceOrientations()[cell_local_id];
    std::vector<double> face_mu_values(cell_num_faces);

    const auto& sigma_t = xs_.at(cell.block_id)->GetSigmaTotal();
    const auto& G = unit_cell_matrices_[cell_local_id].intV_shapeI_gradshapeJ;
    const auto& M = unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;
    const auto& M_surf = unit_cell_matrices_[cell_local_id].intS_shapeI_shapeJ;
    const auto& M_r = secondary_unit_cell_matrices_[cell_local_id].intV_shapeI_shapeJ;

    const auto ni_deploc_face_counter = deploc_face_counter;
    const auto ni_preloc_face_counter = preloc_face_counter;
    for (size_t as_ss_idx = 0; as_ss_idx < as_angle_indices.size(); ++as_ss_idx)
    {
      const auto direction_num = as_angle_indices[as_ss_idx];
      const auto& omega = groupset_.quadrature->GetOmega(direction_num);
      const auto wt = groupset_.quadrature->GetWeight(direction_num);
      const double tau = fac_diamond_difference[direction_num];
      const double c_ang = fac_streaming_operator[direction_num];
      // The zero-weight final direction mu = +1 ends the recursion. Its flux is the angular-edge
      // flux psi_{N+1/2}; solving for it would be singular in a void cell at the center, whose
      // inflow face has zero area. It contributes nothing to the flux moments and is used only as
      // reflected inflow for the starting direction.
      const bool final_direction = wt == 0.0 and omega.z > 0.0;

      deploc_face_counter = ni_deploc_face_counter;
      preloc_face_counter = ni_preloc_face_counter;

      // Angular-redistribution term: c_ang (psi_n - psi_{n-1/2}) / r
      for (size_t gsg = 0; gsg < gs_size; ++gsg)
        b[gsg] = Vector<double>(cell_num_nodes, 0.0);
      for (size_t i = 0; i < cell_num_nodes; ++i)
        for (size_t j = 0; j < cell_num_nodes; ++j)
        {
          const auto jr = discretization_.MapDOFLocal(cell, j, unknown_manager_, 0, 0);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            b[gsg](i) += c_ang * M_r(i, j) * psi_edge_[jr + gsg];
        }

      for (size_t i = 0; i < cell_num_nodes; ++i)
        for (size_t j = 0; j < cell_num_nodes; ++j)
          Amat(i, j) = omega.Dot(G(i, j)) + c_ang * M_r(i, j);

      for (size_t f = 0; f < cell_num_faces; ++f)
        face_mu_values[f] = omega.Dot(cell.faces[f].normal);

      // Incoming surface terms
      int in_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::INCOMING)
          continue;

        const auto& cell_face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not cell_face.has_neighbor;
        // The center of a solid sphere has zero area and no incident flux.
        const bool is_center = is_boundary_face and std::fabs(cell_face.centroid.z) < 1.0e-12;

        if (is_local_face)
          ++in_face_counter;
        else if (not is_boundary_face)
          ++preloc_face_counter;

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int j = cell_mapping.MapFaceNode(f, fj);
            const double mu_Nij = -face_mu_values[f] * M_surf[f](i, j);
            Amat(i, j) += mu_Nij;

            const double* psi = nullptr;
            if (is_local_face)
              psi = fluds.UpwindPsi(spls_index, in_face_counter, fj, 0, as_ss_idx);
            else if (not is_boundary_face)
              psi = fluds.NLUpwindPsi(preloc_face_counter, fj, 0, as_ss_idx);
            else if (not is_center)
              psi = angle_set.PsiBoundary(cell_face.neighbor_id,
                                          direction_num,
                                          cell_local_id,
                                          f,
                                          fj,
                                          0,
                                          IsSurfaceSourceActive());
            if (not psi)
              continue;

            for (size_t gsg = 0; gsg < gs_size; ++gsg)
              b[gsg](i) += psi[gsg] * mu_Nij;
          }
        }
      }

      const auto row_offset = static_cast<size_t>(direction_num) * num_moments_;
      const double* m2d_row = m2d_op.data() + row_offset;
      const double* d2m_row = d2m_op.data() + row_offset;

      for (size_t gsg = 0; gsg < gs_size && final_direction; ++gsg)
        for (size_t i = 0; i < cell_num_nodes; ++i)
          b[gsg](i) = psi_edge_[discretization_.MapDOFLocal(cell, i, unknown_manager_, 0, 0) + gsg];

      for (size_t gsg = 0; gsg < gs_size && not final_direction; ++gsg)
      {
        const double sigma_tg = sigma_t[gs_gi + gsg];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp_src = 0.0;
          for (unsigned int m = 0; m < num_moments_; ++m)
            temp_src += m2d_row[m] * source_moments_[cell_transport_view.MapDOF(i, m, gs_gi + gsg)];
          source[i] = temp_src;
        }
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          double temp = 0.0;
          for (size_t j = 0; j < cell_num_nodes; ++j)
          {
            Atemp(i, j) = Amat(i, j) + M(i, j) * sigma_tg;
            temp += M(i, j) * source[j];
          }
          b[gsg](i) += temp;
        }
        GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes));
      }

      // Flux moments
      for (unsigned int m = 0; m < num_moments_; ++m)
      {
        const double wn_d2m = d2m_row[m];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = cell_transport_view.MapDOF(i, m, gs_gi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }

      if (SaveAngularFluxEnabled())
      {
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const size_t imap =
            i * groupset_angle_group_stride_ + direction_num * groupset_group_stride_;
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_psi_data[imap + gsg] = b[gsg](i);
        }
      }

      // Outgoing faces: downstream fluxes, reflected fluxes, and outflow tallies
      int out_face_counter = -1;
      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        if (face_orientations[f] != FaceOrientation::OUTGOING)
          continue;

        out_face_counter++;
        const auto& face = cell.faces[f];
        const bool is_local_face = cell_transport_view.IsFaceLocal(f);
        const bool is_boundary_face = not face.has_neighbor;
        const bool is_reflecting_boundary_face =
          (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());
        const auto& IntF_shapeI = unit_cell_matrices_[cell_local_id].intS_shapeI[f];

        if (not is_boundary_face and not is_local_face)
          ++deploc_face_counter;

        const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            cell_outflow_view.Add(
              f, gs_gi + gsg, wt * face_mu_values[f] * b[gsg](i) * IntF_shapeI(i));

          double* psi = nullptr;
          if (is_local_face)
            psi = fluds.OutgoingPsi(spls_index, out_face_counter, fi, as_ss_idx);
          else if (not is_boundary_face)
            psi = fluds.NLOutgoingPsi(deploc_face_counter, fi, as_ss_idx);
          else if (is_reflecting_boundary_face)
            psi = angle_set.PsiReflected(face.neighbor_id, direction_num, cell_local_id, f, fi);
          else
            continue;

          for (size_t gsg = 0; gsg < gs_size; ++gsg)
            psi[gsg] = b[gsg](i);
        }
      }

      // Angular-edge flux for the next direction: psi_{n+1/2} = (psi_n - (1 - tau) psi_{n-1/2}) /
      // tau. The starting direction (tau = 1) initializes it.
      const double f0 = 1.0 / tau;
      const double f1 = f0 - 1.0;
      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = discretization_.MapDOFLocal(cell, i, unknown_manager_, 0, 0);
        for (size_t gsg = 0; gsg < gs_size; ++gsg)
          psi_edge_[ir + gsg] = f0 * b[gsg](i) - f1 * psi_edge_[ir + gsg];
      }
    }
  }
}

} // namespace opensn
