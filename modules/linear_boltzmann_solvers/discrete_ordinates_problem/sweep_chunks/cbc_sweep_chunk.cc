// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/cell/cell.h"
#include "framework/logging/log.h"
#include "caliper/cali.h"

namespace opensn
{

/**
 * @brief Constructs a CBCSweepChunk object.
 *
 * Initializes base class and member variables. Strides related to angular flux data
 * are initialized later in `SetAngleSet`.
 *
 * @param destination_phi Reference to the global vector for storing scalar flux moments.
 * @param destination_psi Reference to the global vector for storing angular fluxes
 *                        (written to if `save_angular_flux` is enabled).
 * @param grid Shared pointer to the mesh continuum.
 * @param discretization Reference to the spatial discretization scheme.
 * @param unit_cell_matrices Vector of precomputed unit cell matrices.
 * @param cell_transport_views Vector of LBS cell views for transport properties.
 * @param densities Vector of cell material densities.
 * @param source_moments Vector of source moments.
 * @param groupset Reference to the LBS groupset being solved.
 * @param xs Map of material IDs to multi-group cross sections.
 * @param num_moments Number of flux moments.
 * @param max_num_cell_dofs Maximum number of degrees of freedom on any cell.
 */
CBCSweepChunk::CBCSweepChunk(std::vector<double>& destination_phi,
                             std::vector<double>& destination_psi,
                             const std::shared_ptr<MeshContinuum> grid,
                             const SpatialDiscretization& discretization,
                             const std::vector<UnitCellMatrices>& unit_cell_matrices,
                             std::vector<CellLBSView>& cell_transport_views,
                             const std::vector<double>& densities,
                             const std::vector<double>& source_moments,
                             const LBSGroupset& groupset,
                             const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                             int num_moments,
                             int max_num_cell_dofs)
  : SweepChunk(destination_phi,
               destination_psi,
               grid,
               discretization,
               unit_cell_matrices,
               cell_transport_views,
               densities,
               source_moments,
               groupset,
               xs,
               num_moments,
               max_num_cell_dofs),
    fluds_(nullptr),
    gs_size_(0),
    gs_gi_(0),
    // Initialize new stride members
    num_angles_in_set_local_(0),
    local_compact_angle_stride_(0),
    local_compact_node_stride_(0),
    num_angles_in_set_remote_(0),
    remote_angle_stride_(0),
    remote_node_stride_(0),
    surface_source_active_(false),
    cell_(nullptr),
    cell_local_id_(0),
    cell_mapping_(nullptr),
    cell_transport_view_(nullptr),
    cell_num_faces_(0),
    cell_num_nodes_(0)
{
}

/**
 * @brief Sets the context for the current `AngleSet` being processed.
 *
 * For each new `AngleSet`, this method:
 * - Stores a pointer to the `CBC_FLUDS` instance associated with the `angle_set`.
 * - Caches the number of groups (`gs_size_`) and the starting global group index (`gs_gi_`)
 *   from the parent `LBSGroupset`.
 * - Determines if surface sources on boundaries are active.
 * - Calculates and stores strides necessary for indexing into:
 *   - The compact `local_psi_data_` within `CBC_FLUDS` (for local cell-to-cell transfers).
 *   - MPI communication buffers (for remote data exchange).
 *
 * @param angle_set Reference to the `AngleSet` that this chunk will operate on.
 */
void
CBCSweepChunk::SetAngleSet(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CbcSweepChunk::SetAngleSet");

  fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());

  /// Number of groups in LBSGroupset
  gs_size_ = groupset_.groups.size();
  gs_gi_ = groupset_.groups.front().id;

  surface_source_active_ = IsSurfaceSourceActive();

  /// --- Strides for `local_psi_data_` (used for local cell-to-cell transfers) ---

  /// Angles specific to this AngleSet
  num_angles_in_set_local_ = angle_set.GetNumAngles();

  /// Groups are contiguous per angle within a spatial_dof-angle block
  local_compact_angle_stride_ = gs_size_;

  /// Stride between spatial_dof blocks
  local_compact_node_stride_ = num_angles_in_set_local_ * gs_size_;

  /// --- Strides for REMOTE communication buffer layout (for psi_dnwnd_data) ---

  /// This buffer is structured: face spatial DOF major -> angle in set major -> group major

  /// Same as local, for clarity of use
  num_angles_in_set_remote_ = angle_set.GetNumAngles();

  /// Stride between angles within a face-node's data block
  remote_angle_stride_ = gs_size_;

  /// Stride between face-node data blocks
  remote_node_stride_ = num_angles_in_set_remote_ * gs_size_;
}

/**
 * @brief Sets the current cell for which the sweep operation will be performed.
 *
 * Caches cell-specific information like its local ID, cell mapping, LBS view,
 * number of faces/nodes, and precomputed unit cell matrices to optimize
 * access during the `Sweep` method.
 *
 * @param cell_ptr Pointer to the `Cell` to be processed.
 * @param angle_set Reference to the current `AngleSet`
 */
void
CBCSweepChunk::SetCell(const Cell* cell_ptr, AngleSet& angle_set)
{
  cell_ = cell_ptr;
  cell_local_id_ = cell_ptr->local_id;
  cell_mapping_ = &discretization_.GetCellMapping(*cell_);
  cell_transport_view_ = &cell_transport_views_[cell_->local_id];
  cell_num_faces_ = cell_->faces.size();
  cell_num_nodes_ = cell_mapping_->GetNumNodes();

  /// Get cell matrices
  G_ = unit_cell_matrices_[cell_local_id_].intV_shapeI_gradshapeJ;
  M_ = unit_cell_matrices_[cell_local_id_].intV_shapeI_shapeJ;
  M_surf_ = unit_cell_matrices_[cell_local_id_].intS_shapeI_shapeJ;
  IntS_shapeI_ = unit_cell_matrices_[cell_local_id_].intS_shapeI;
}

/**
 * @brief Executes the discrete ordinates sweep calculation for the current cell.
 *
 * This method iterates over all angles in the current `AngleSet` (via `as_ss_idx`)
 * and all energy groups in the `LBSGroupset` (via `gsg`). For each angle/group:
 * 1. Assembles the local linear system representing the transport equation
 *    (streaming, collision, surface terms).
 * 2. Incorporates upwind angular fluxes from:
 *    - Local neighbor cells: Reads from `fluds_->local_psi_data_` using compact strides
 *      and `as_ss_idx`.
 *    - Remote (MPI) neighbor cells: Reads from MPI data packets managed by `fluds_`,
 *      interpreting them using `as_ss_idx` and remote strides.
 *    - Boundaries: Uses `angle_set.PsiBoundary()`.
 * 3. Adds volumetric source contributions (fixed, scattering, fission).
 * 4. Solves the local system for the cell's nodal angular fluxes.
 * 5. Updates the global scalar flux moments (`destination_phi_`).
 * 6. If `save_angular_flux_` is true, saves the computed angular fluxes to the
 *    global `destination_psi_` vector (using global angle indices and full strides).
 * 7. Propagates outgoing angular fluxes:
 *    - To local downwind cells: Writes to `fluds_->local_psi_data_` using compact
 *      strides and `as_ss_idx`.
 *    - To remote downwind cells: Stages data into MPI send buffers via `async_comm_`.
 *    - To reflecting boundaries: Uses `angle_set.PsiReflected()`.
 *    - Updates outflow tallies for particle balance.
 *
 * @param angle_set Reference to the current `AngleSet` being processed.
 */
void
CBCSweepChunk::Sweep(AngleSet& angle_set)
{
  const auto& m2d_op = groupset_.quadrature->GetMomentToDiscreteOperator();
  const auto& d2m_op = groupset_.quadrature->GetDiscreteToMomentOperator();

  /// Initialize local matrices and vectors (sized appropriately for the current cell)
  DenseMatrix<double> Amat(max_num_cell_dofs_, max_num_cell_dofs_);
  DenseMatrix<double> Atemp(max_num_cell_dofs_, max_num_cell_dofs_);
  std::vector<Vector<double>> b(gs_size_, Vector<double>(max_num_cell_dofs_));
  std::vector<double> source(max_num_cell_dofs_);

  const auto& face_orientations = angle_set.GetSPDS().GetCellFaceOrientations()[cell_local_id_];
  std::vector<double> face_mu_values(cell_num_faces_);

  const auto& rho = densities_[cell_local_id_];
  const auto& sigma_t = xs_.at(cell_->block_id)->GetSigmaTotal();

  /// Global angle indices for this set
  const std::vector<size_t>& as_angle_indices = angle_set.GetAngleIndices();

  /// Loop over angles in *THIS* AngleSet using a LOCAL index (as_ss_idx)
  for (size_t as_ss_idx = 0; as_ss_idx < num_angles_in_set_local_; ++as_ss_idx)
  {
    /// Get the GLOBAL angle index
    auto direction_num = as_angle_indices[as_ss_idx];
    auto omega = groupset_.quadrature->omegas[direction_num];
    auto wt = groupset_.quadrature->weights[direction_num];

    /// Reset RHS vector b for the current angle
    for (int gsg = 0; gsg < gs_size_; ++gsg)
      for (int i = 0; i < cell_num_nodes_; ++i)
        b[gsg](i) = 0.0;

    /// Assemble streaming term Amat = omega dot Grad(Shape)
    for (int i = 0; i < cell_num_nodes_; ++i)
      for (int j = 0; j < cell_num_nodes_; ++j)
        Amat(i, j) = omega.Dot(G_(i, j));

    /// Calculate mu_n for each face for the current angle
    for (int f = 0; f < cell_num_faces_; ++f)
      face_mu_values[f] = omega.Dot(cell_->faces[f].normal);

    /// Surface integrals (Upstream contributions to Amat and RHS b)
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::INCOMING)
        continue;

      const auto& face = cell_->faces[f];
      const bool is_local_face = cell_transport_view_->IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      auto face_nodal_mapping = &fluds_->GetCommonData().GetFaceNodalMapping(cell_local_id_, f);

      /// For remote faces, get the pre-received data block (contains all angles in set for this
      /// face)
      const std::vector<double>* psi_nonlocal_upwnd_data_block = nullptr;
      if ((not is_local_face) and (not is_boundary_face))
      {
        psi_nonlocal_upwnd_data_block = &fluds_->GetNonLocalUpwindData(cell_->global_id, f);
      }

      const size_t num_face_nodes = cell_mapping_->GetNumFaceNodes(f);

      /// Loop over nodes ON THIS FACE of current cell
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        /// Map face node fi to cell node i
        const int i = cell_mapping_->MapFaceNode(f, fi);

        /// Loop over basis functions on this face
        for (int fj = 0; fj < num_face_nodes; ++fj)
        {
          /// Map face node fj to cell node j
          const int j = cell_mapping_->MapFaceNode(f, fj);

          const double mu_Nij = -face_mu_values[f] * M_surf_[f](i, j);
          Amat(i, j) += mu_Nij; /// Add to LHS matrix

          const double* psi_upwind_groups_ptr = nullptr;

          if (is_local_face)
          {
            /// Local upwind read
            const Cell* upwind_cell = cell_transport_view_->FaceNeighbor(f);

            /// Node index WITHIN the UPWIND cell corresponding to face node fj of CURRENT cell's
            /// face f
            const unsigned int adj_cell_node_idx = face_nodal_mapping->cell_node_mapping_[fj];

            /// 1. Get base pointer to the start of upwind_cell's data block in compact
            /// local_psi_data_
            const double* psi_upwind_cell_base_ptr = fluds_->GetLocalUpwindPsi(*upwind_cell);

            /// 2. Calculate relative offset to the specific node and angle within that block
            const size_t offset_in_cell_block =
              adj_cell_node_idx *
                local_compact_node_stride_ + /// Offset to this node's (all angles) block
              as_ss_idx *
                local_compact_angle_stride_; /// Further offset to this specific angle's group block

            psi_upwind_groups_ptr = psi_upwind_cell_base_ptr + offset_in_cell_block;
          }
          else if (not is_boundary_face) /// Remote face
          {
            /// Remote upwind read (interprets multi-AngleSet-angle packet)
            assert(psi_nonlocal_upwnd_data_block != nullptr);
            const unsigned int adj_face_node_idx = face_nodal_mapping->face_node_mapping_[fj];
            psi_upwind_groups_ptr = fluds_->GetNonLocalUpwindPsi(
              *psi_nonlocal_upwnd_data_block, adj_face_node_idx, as_ss_idx);
          }
          else /// Physical boundary face
          {
            /// Boundary upwind read
            psi_upwind_groups_ptr = angle_set.PsiBoundary(face.neighbor_id,
                                                          direction_num,
                                                          cell_local_id_,
                                                          f,
                                                          fj,
                                                          gs_gi_,
                                                          surface_source_active_);
          }

          if (psi_upwind_groups_ptr != nullptr)
          {
            for (int gsg = 0; gsg < gs_size_; ++gsg)
            {
              b[gsg](i) += psi_upwind_groups_ptr[gsg] * mu_Nij;
            }
          }
        } /// for fj (face basis function)
      }   /// for fi (face node)
    }     /// for f (face)

    /// Volumetric source term (Q) and Mass term (Sigma_t * M * psi) assembly
    for (int gsg = 0; gsg < gs_size_; ++gsg)
    {
      /// Total XS for global group gs_gi_ + gsg
      double sigma_tg = rho * sigma_t[gs_gi_ + gsg];

      /// Assemble nodal source `source[i]` from `source_moments_` for the current angle
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp_src = 0.0;
        for (int m = 0; m < num_moments_; ++m)
        {
          const size_t ir = cell_transport_view_->MapDOF(i, m, static_cast<int>(gs_gi_ + gsg));
          temp_src += m2d_op[m][direction_num] * source_moments_[ir];
        }
        source[i] = temp_src; /// Source for node i, current angle, group gsg
      }

      /// Add Mass matrix term to Atemp and source term to b
      /// Atemp = Amat + sigma_tgr * M
      /// b += M * Q
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        double temp_M_source = 0.0;
        for (int j = 0; j < cell_num_nodes_; ++j)
        {
          const double Mij = M_(i, j); /// Mass matrix element
          Atemp(i, j) = Amat(i, j) + Mij * sigma_tg;
          temp_M_source += Mij * source[j];
        }
        b[gsg](i) += temp_M_source; /// Add volumetric source contribution to RHS
      }
      /// Solve the local system Atemp * psi_nodal_gsg = b_gsg for psi_nodal_gsg
      /// (solution stored in b[gsg])
      GaussElimination(Atemp, b[gsg], static_cast<int>(cell_num_nodes_));
    } /// for gsg (group in set)

    /// Update scalar flux moments (phi) using the solved angular flux (b[gsg])
    for (int m = 0; m < num_moments_; ++m)
    {
      /// Quadrature weight * D2M operator
      const double wn_d2m = d2m_op[m][direction_num];
      for (int i = 0; i < cell_num_nodes_; ++i)
      {
        /// Map to global phi vector
        const size_t ir = cell_transport_view_->MapDOF(i, m, gs_gi_);
        for (int gsg = 0; gsg < gs_size_; ++gsg)
        {
          /// b[gsg](i) is solved angular flux for node i, group gsg
          destination_phi_[ir + gsg] += wn_d2m * b[gsg](i);
        }
      }
    }

    /// If requested, save angular fluxes during sweep
    /// Crucically, this *does not* write to local_psi_data_,
    /// but instead to `psi_new_local_[groupset.id]` via reference
    if (save_angular_flux_)
    {
      // Get a raw pointer to the beginning of the current cell's data block within
      // destination_psi_
      // discretization_.MapDOFLocal maps the cell's node 0, using the LBSGroupset's
      // full psi_uk_man_ (which understands all angles in the groupset's
      // quadrature), to get the starting 1D index.
      double* cell_psi_data_base_ptr =
        &destination_psi_[discretization_.MapDOFLocal(*cell_, 0, groupset_.psi_uk_man_, 0, 0)];

      /// Iterate over each node of the current cell
      for (size_t i = 0; i < cell_num_nodes_; ++i)
      {
        /// Calculate the offset within the cell's angular flux data block to reach the
        /// data for the current node 'i' and the current GLOBAL angle 'direction_num'
        ///
        /// - groupset_angle_group_stride_:
        ///   Stride to get from one node's full data block (all global angles,
        ///   all groups) to the next node's block
        ///   This stride is based on the total number of angles in the
        ///   groupset's quadrature
        ///
        /// - groupset_group_stride_:
        ///   Stride to get from one global angle's group data to the next global
        ///   angle's group data (for the same node)
        ///   This is gs_size_
        const size_t relative_offset_in_cell =
          i * groupset_angle_group_stride_ + /// Offset to current node 'i'
          direction_num *
            groupset_group_stride_; /// Further offset to current GLOBAL angle 'direction_num'

        /// Iterate over each group in the current groupset
        for (int gsg = 0; gsg < gs_size_; ++gsg)
        {
          /// Write the solved angular flux (b[gsg](i) for current cell node 'i',
          /// current angle, and group 'gsg') into the correct location in the
          /// output_psi vector
          /// The final index is base_ptr + relative_offset + group_index
          cell_psi_data_base_ptr[relative_offset_in_cell + gsg] = b[gsg](i);
        }
      }
    }

    /// Perform outgoing surface operations
    /// (write to local/remote/reflecting boundaries)
    for (int f = 0; f < cell_num_faces_; ++f)
    {
      if (face_orientations[f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face = cell_->faces[f];
      const size_t num_face_nodes = cell_mapping_->GetNumFaceNodes(f);
      const bool is_local_face = cell_transport_view_->IsFaceLocal(f);
      const bool is_boundary_face = not face.has_neighbor;
      const bool is_reflecting_boundary_face =
        (is_boundary_face and angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting());

      /// Surface integral of shape functions for face f
      const auto& IntFi_shapeI_vec = IntS_shapeI_[f];

      /// Prepare buffer for remote downwind data
      std::vector<double>* psi_nonlocal_dnwnd_data_block_ptr = nullptr;

      if (not is_boundary_face and not is_local_face)
      { /// If remote face
        const int locality = cell_transport_view_->FaceLocality(f);
        auto& face_nodal_mapping = fluds_->GetCommonData().GetFaceNodalMapping(cell_local_id_, f);
        auto& async_comm = *angle_set.GetCommunicator();

        /// Remote downwind buffer sizing (for multi-AngleSet-angle packet)
        /// Size of buffer for ALL angles in this AngleSet, ALL groups,
        /// for ALL nodes on this face
        /// remote_node_stride_ = num_angles_in_set_remote_ * gs_size_
        const size_t data_size_for_msg = num_face_nodes * remote_node_stride_;

        psi_nonlocal_dnwnd_data_block_ptr = &async_comm.InitGetDownwindMessageData(
          locality,
          face.neighbor_id,
          face_nodal_mapping.associated_face_,
          angle_set.GetID(),
          data_size_for_msg); /// Request buffer for multi-angle packet
      }

      /// Loop over nodes ON THIS FACE of current cell
      for (int fi = 0; fi < num_face_nodes; ++fi)
      {
        /// Node index in current cell
        const int i = cell_mapping_->MapFaceNode(f, fi);

        /// Tally outflow for particle balance (for ALL boundary faces)
        if (is_boundary_face)
        {
          /// Integral of shape func for node i on face f
          const double IntFi_shapeI_val = IntFi_shapeI_vec(i);
          for (int gsg = 0; gsg < gs_size_; ++gsg)
            cell_transport_view_->AddOutflow(
              f, gs_gi_ + gsg, wt * face_mu_values[f] * b[gsg](i) * IntFi_shapeI_val);
        }

        /// Pointer to start of group data for outgoing psi
        double* psi_downwind_groups_ptr = nullptr;

        if (is_local_face)
        {
          /// Local downwind write
          /// 1. Get base pointer to the start of current_cell's data block in
          /// compact local_psi_data_
          double* psi_downwind_cell_base_ptr = fluds_->GetLocalDownwindPsi(*cell_);

          /// 2. Calculate relative offset to the specific node and angle within
          /// that block
          const size_t offset_in_cell_block =
            i * local_compact_node_stride_ + /// Offset to this node's (all angles) block
            as_ss_idx *
              local_compact_angle_stride_; /// Further offset to this specific angle's group block

          psi_downwind_groups_ptr = psi_downwind_cell_base_ptr + offset_in_cell_block;
        }
        else if (not is_boundary_face) /// Remote face
        {
          /// Remote downwind write (indexes into multi-AngleSet-angle packet)
          assert(psi_nonlocal_dnwnd_data_block_ptr != nullptr);

          /// Indexing into the psi_dnwnd_data buffer
          /// This buffer has layout: spatial DOF major -> angle in set major -> group major
          /// fi = current face node index (0 to num_face_nodes_on_this_face-1)
          /// as_ss_idx = current AngleSet angle index (0 to num_angles_in_set_remote_-1)
          const size_t nonlocal_addr_offset =
            fi * remote_node_stride_ + /// Offset to the block for this face node
            as_ss_idx *
              remote_angle_stride_; /// Offset to the block for this angle within the node's block

          psi_downwind_groups_ptr = &(*psi_nonlocal_dnwnd_data_block_ptr)[nonlocal_addr_offset];
        }
        else if (is_reflecting_boundary_face)
        {
          psi_downwind_groups_ptr =
            angle_set.PsiReflected(face.neighbor_id,
                                   direction_num, /// direction_num is global angle index
                                   cell_local_id_,
                                   f,
                                   fi);
        }

        /// Write the solved angular flux (b[gsg](i)) to the determined location
        if (psi_downwind_groups_ptr != nullptr)
        {
          for (int gsg = 0; gsg < gs_size_; ++gsg)
          {
            psi_downwind_groups_ptr[gsg] = b[gsg](i);
          }
        }
      } /// for fi (face node)
    }   /// for f (outgoing face)
  }     /// for as_ss_idx (angle in set)
}

} // namespace opensn
