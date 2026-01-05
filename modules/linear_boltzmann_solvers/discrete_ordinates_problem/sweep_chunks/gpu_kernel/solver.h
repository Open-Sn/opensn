// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/mesh_view.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/quadrature_view.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn::gpu_kernel
{

/// Compute the sweep matrix from gradient, mass and the source term
template <std::size_t ndofs>
__device__ void
ComputeGMS(double* sweep_matrix,
           double* psi,
           double* s,
           CellView& cell,
           DirectionView& direction,
           const unsigned int& group_idx,
           const std::uint32_t& num_moments,
           const Arguments& args)
{
  // get sigmaT
  double sigma_t = cell.total_xs[args.groupset_start + group_idx];
  // compute source term
  const double* src_moment = args.src_moment + cell.phi_address + args.groupset_start + group_idx;
  _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
  {
    double src_per_moment = 0.0;
    for (std::uint32_t m = 0; m < num_moments; ++m)
    {
      src_per_moment += direction.m2d[m] * (*src_moment);
      src_moment += args.num_groups;
    }
    s[i] = src_per_moment;
  }
  // add source, transfer and mass contribution
  double* A = sweep_matrix;
  const std::array<double, 4>* GM_data =
    reinterpret_cast<const std::array<double, 4>*>(cell.GM_data);
  _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
  {
    _Pragma("unroll") for (std::uint32_t j = 0; j < ndofs; ++j)
    {
      std::array<double, 4> GM = *(GM_data++);
      // compute A += G * Omega + M * sigma_t
      A[j] += direction.omega[0] * GM[0] + direction.omega[1] * GM[1] + direction.omega[2] * GM[2] +
              sigma_t * GM[3];
      // compute psi += M @ s
      psi[i] += GM[3] * s[j];
    }
    A += ndofs;
  }
}

/// Compute the sweep matrix from surface integral
template <std::size_t ndofs>
__device__ void
ComputeSurfaceIntegral(double* sweep_matrix,
                       double* psi,
                       CellView& cell,
                       DirectionView& direction,
                       const std::uint64_t* cell_edge_data,
                       const unsigned int& angle_group_idx,
                       const Arguments& args)
{
  // loop over each face
  std::uint32_t face_node_counter = 0;
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    // get face view
    FaceView face;
    cell.GetFaceView(face, f);
    // determine if this face is incoming
    AAHD_NodeIndex idx(cell_edge_data[face_node_counter]);
    if (idx.IsUndefined() || idx.IsOutGoing())
    {
      face_node_counter += face.num_face_nodes;
      continue;
    }
    double mu = direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
                direction.omega[2] * face.normal[2];
    // compute surface integral
    for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
    {
      std::uint32_t i = face.cell_mapping_data[fi];
      double* Ai = sweep_matrix + i * ndofs;
      for (std::uint32_t fj = 0; fj < face.num_face_nodes; ++fj)
      {
        std::uint32_t j = face.cell_mapping_data[fj];
        double mu_Nij = -mu * face.M_surf_data[fi * face.num_face_nodes + fj];
        Ai[j] += mu_Nij;
        double* upwind_psi =
          args.flud_data.GetIncomingFluxPointer(cell_edge_data[face_node_counter + fj]);
        psi[i] += upwind_psi[angle_group_idx] * mu_Nij;
      }
    }
    // update face node counter
    face_node_counter += face.num_face_nodes;
  }
}

/// Gaussian elimination
template <std::size_t ndofs>
__device__ void
GaussianElimination(double* sweep_matrix, double* psi)
{
  // forward elimination
  double* A_i = sweep_matrix;
  _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
  {
    double inv_diag = 1.0 / A_i[i];
    // normalize the pivot row
    _Pragma("unroll") for (std::uint32_t j = i; j < ndofs; ++j)
    {
      A_i[j] *= inv_diag;
    }
    psi[i] *= inv_diag;
    // eliminate rows below
    double* A_k = A_i + ndofs;
    _Pragma("unroll") for (std::uint32_t k = i + 1; k < ndofs; ++k)
    {
      double factor = -A_k[i];
      _Pragma("unroll") for (std::uint32_t j = i; j < ndofs; ++j)
      {
        A_k[j] += factor * A_i[j];
      }
      psi[k] += factor * psi[i];
      A_k += ndofs;
    }
    A_i += ndofs;
  }
  // back substitution — row-wise access
  if constexpr (ndofs >= 2)
  {
    _Pragma("unroll") for (std::int32_t j = ndofs - 2; j >= 0; --j)
    {
      double* A_j = sweep_matrix + j * ndofs;
      _Pragma("unroll") for (std::int32_t i = j + 1; i < ndofs; ++i)
      {
        psi[j] -= A_j[i] * psi[i];
      }
    }
  }
}

/// Record angular flux to downwind and compute outflow for boundary faces.
__device__ inline void
WritePsiToFludsAndOutflow(double* psi,
                          CellView& cell,
                          DirectionView& direction,
                          const std::uint64_t* cell_edge_data,
                          const unsigned int& angle_group_idx,
                          const unsigned int& group_idx,
                          const Arguments& args)
{
  // loop over each face
  std::uint32_t face_node_counter = 0;
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    // get face view
    FaceView face;
    cell.GetFaceView(face, f);
    // determine if this face is outgoing
    AAHD_NodeIndex idx(cell_edge_data[face_node_counter]);
    if (idx.IsUndefined() || !idx.IsOutGoing())
    {
      face_node_counter += face.num_face_nodes;
      continue;
    }
    double mu = direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
                direction.omega[2] * face.normal[2];
    // loop over each face node
    for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
    {
      std::uint32_t i = face.cell_mapping_data[fi];
      // put copy psi to FLUD
      double* downwind_psi =
        args.flud_data.GetOutgoingFluxPointer(cell_edge_data[face_node_counter + fi]);
      downwind_psi[angle_group_idx] = psi[i];
      // compute ouflow for boundary face
      if (face.outflow != nullptr)
      {
        double outflow = direction.weight * mu * face.IntS_shapeI_data[fi] * psi[i];
        crb::atomic_add(face.outflow + args.groupset_start + group_idx, outflow);
      }
    }
    face_node_counter += face.num_face_nodes;
  }
}

/// Compute the scalar flux.
template <std::size_t ndofs>
__device__ void
ComputePhi(double* psi,
           CellView& cell,
           DirectionView& direction,
           const unsigned int& group_idx,
           const std::uint32_t& num_moments,
           const Arguments& args)
{
  double* phi = args.phi + cell.phi_address + args.groupset_start + group_idx;
  _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
  {
    for (std::uint32_t m = 0; m < num_moments; ++m)
    {
      crb::atomic_add(phi, direction.d2m[m] * psi[i]);
      phi += args.num_groups;
    }
  }
}

/// Store the angular flux
template <std::size_t ndofs>
__device__ void
SaveAngularFlux(const double* psi,
                double* saved_psi,
                CellView& cell,
                const unsigned int& angle_group_idx,
                const std::uint64_t& stride_size)
{
  saved_psi += cell.save_psi_index * stride_size + angle_group_idx;
  _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
  {
    *saved_psi = psi[i];
    saved_psi += stride_size;
  }
}

/// Template device function performing the sweep
template <std::size_t ndofs>
__device__ void
Sweep(const Arguments& args,
      CellView& cell,
      DirectionView& direction,
      const std::uint64_t* cell_edge_data,
      const unsigned int& angle_group_idx,
      const unsigned int& group_idx,
      const std::uint32_t& num_moments,
      double* saved_psi)
{
  // initialize buffer
  Buffer<ndofs> buffer;
  // prepare linear system to solve
  ComputeGMS<ndofs>(
    buffer.A(), buffer.b(), buffer.s(), cell, direction, group_idx, num_moments, args);
  ComputeSurfaceIntegral<ndofs>(
    buffer.A(), buffer.b(), cell, direction, cell_edge_data, angle_group_idx, args);
  // solve for the angular flux
  GaussianElimination<ndofs>(buffer.A(), buffer.b());
  // save the result
  WritePsiToFludsAndOutflow(
    buffer.b(), cell, direction, cell_edge_data, angle_group_idx, group_idx, args);
  ComputePhi<ndofs>(buffer.b(), cell, direction, group_idx, num_moments, args);
  if (saved_psi != nullptr)
  {
    SaveAngularFlux<ndofs>(
      buffer.b(), saved_psi, cell, angle_group_idx, args.flud_data.stride_size);
  }
}

} // namespace opensn::gpu_kernel
