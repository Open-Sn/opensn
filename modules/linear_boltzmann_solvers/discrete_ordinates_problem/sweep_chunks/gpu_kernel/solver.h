// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/buffer.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/mesh_view.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/quadrature_view.h"
#include "caribou/main.hpp"
#include <array>
#include <cstdint>
#include <type_traits>

namespace crb = caribou;

namespace opensn::gpu_kernel
{

// NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast, modernize-use-auto)

/**
 * Default maximum DOF count for register-resident local systems.
 */
constexpr std::size_t max_register_resident_dofs = 4;

/**
 * Add a value to an array entry whose index is known only at run time.
 * For systems of up to max_register_resident_dofs DOFs, every access uses a compile-time index.
 * The local array can then stay in registers instead of being placed in local memory.
 *
 * \tparam n Array length.
 * \tparam ndofs Cell DOF count.
 * \param array Array to update.
 * \param index Index of the entry.
 * \param value Value to add.
 */
template <std::size_t n, std::size_t ndofs = n>
__CRB_DEVICE_FUNC__ void
AddAtIndex(double* array, std::uint32_t index, double value)
{
  if constexpr (ndofs <= max_register_resident_dofs)
  {
    _Pragma("unroll") for (std::uint32_t k = 0; k < n; ++k)
    {
      if (k == index)
        array[k] += value;
    }
  }
  else
    array[index] += value;
}

/**
 * Return an array entry whose index is known only at run time (see AddAtIndex).
 *
 * \tparam n Array length.
 * \tparam register_resident Whether the array is kept in registers (compile-time indices only).
 * \tparam T Entry type.
 * \param array Array to read.
 * \param index Index of the entry.
 * \return Value of the entry.
 */
template <std::size_t n, bool register_resident = (n <= max_register_resident_dofs), typename T>
__CRB_DEVICE_FUNC__ T
ValueAtIndex(const T* array, std::uint32_t index)
{
  if constexpr (register_resident)
  {
    T value{};
    _Pragma("unroll") for (std::uint32_t k = 0; k < n; ++k)
    {
      if (k == index)
        value = array[k];
    }
    return value;
  }
  else
    return array[index];
}

/// Compute the sweep matrix from gradient, mass, and the source term
template <std::size_t ndofs, SweepKind k>
__CRB_DEVICE_FUNC__ void
ComputeGMS(double* sweep_matrix,
           double* psi,
           double* s,
           CellView& cell,
           DirectionView& direction,
           const unsigned int& group_idx,
           const std::uint32_t& num_moments,
           const Arguments<k>& args)
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

/// Compute the direction cosine. On AMD GPUs, round the y product separately to match
/// run-time-indexed assembly.
__CRB_DEVICE_FUNC__ double
CosineRoundingSecondProduct(const DirectionView& direction, const FaceView& face)
{
#if defined(__HIP_DEVICE_COMPILE__) && defined(__AMDGCN__)
  return __builtin_fma(
    direction.omega[2],
    face.normal[2],
    __builtin_fma(direction.omega[0], face.normal[0], direction.omega[1] * face.normal[1]));
#else
  return direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
         direction.omega[2] * face.normal[2];
#endif
}

/// Add without contracting the sum with a product used as an addend.
__CRB_DEVICE_FUNC__ double
AddRounded(double a, double b)
{
#ifdef __CUDA_ARCH__
  return __dadd_rn(a, b);
#else
#if defined(__clang__) && !defined(__NVCC__)
#pragma clang fp contract(off)
#endif
  return a + b;
#endif
}

/// Match run-time-indexed matrix assembly. AMD GPUs fuse -mu * N with the sum.
/// Other backends add the separately rounded product mu_Nij.
__CRB_DEVICE_FUNC__ double
AddFaceCouplingToEntry(double a, double mu, double N, double mu_Nij)
{
#if defined(__HIP_DEVICE_COMPILE__) && defined(__AMDGCN__)
  static_cast<void>(mu_Nij);
  return __builtin_fma(N, -mu, a);
#else
  static_cast<void>(mu);
  static_cast<void>(N);
  return AddRounded(a, mu_Nij);
#endif
}

/// Assemble an incoming face for a register-resident CBC system above the default DOF limit.
/// Face nodes must map to distinct cell nodes, so each matrix entry has at most one contribution
/// per face. Apply right-hand-side contributions in face-node order.
template <std::size_t ndofs, bool register_resident = (ndofs <= max_register_resident_dofs)>
__CRB_DEVICE_FUNC__ void
AddIncomingFaceCBC(double* sweep_matrix,
                   double* psi,
                   const FaceView& face,
                   const double mu,
                   const std::uint64_t* face_edge_data,
                   const unsigned int& angle_group_idx,
                   const unsigned int& group_idx,
                   const Arguments<SweepKind::CBC>& args)
{
  static_assert(register_resident and ndofs > max_register_resident_dofs);
  const std::uint32_t num_face_nodes = face.num_face_nodes;
  // The node indices are loaded before any is used, so that their loads overlap.
  std::array<std::uint32_t, ndofs> node{};
  std::array<std::uint64_t, ndofs> node_index{};
  _Pragma("unroll") for (std::uint32_t fj = 0; fj < ndofs; ++fj)
  {
    if (fj < num_face_nodes)
    {
      node[fj] = face.cell_mapping_data[fj];
      node_index[fj] = face_edge_data[fj];
    }
  }
  std::array<double, ndofs> upwind_psi{};
  _Pragma("unroll") for (std::uint32_t fj = 0; fj < ndofs; ++fj)
  {
    if (fj < num_face_nodes)
      upwind_psi[fj] = *args.flud_data.GetIncomingFluxPointer(
        node_index[fj], angle_group_idx, group_idx, args.boundary, args.is_surface_source_active);
  }
  if constexpr (register_resident and ndofs > max_register_resident_dofs)
  {
    // Find the face node from its compile-time column so the matrix stays in registers.
    std::array<std::uint32_t, ndofs> face_node_of{};
    _Pragma("unroll") for (std::uint32_t c = 0; c < ndofs; ++c)
    {
      face_node_of[c] = num_face_nodes;
      _Pragma("unroll") for (std::uint32_t fj = 0; fj < ndofs; ++fj)
      {
        if (fj < num_face_nodes and node[fj] == c)
          face_node_of[c] = fj;
      }
    }
    _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
    {
      const std::uint32_t fi = face_node_of[i];
      if (fi == num_face_nodes)
        continue;
      const double* M_surf_row = face.M_surf_data + static_cast<std::size_t>(fi * num_face_nodes);
      std::array<double, ndofs> M_surf{};
      _Pragma("unroll") for (std::uint32_t fj = 0; fj < ndofs; ++fj)
      {
        if (fj < num_face_nodes)
          M_surf[fj] = M_surf_row[fj];
      }
      double* Ai = sweep_matrix + i * ndofs;
      _Pragma("unroll") for (std::uint32_t c = 0; c < ndofs; ++c)
      {
        const std::uint32_t fj = face_node_of[c];
        double N = M_surf[0];
        _Pragma("unroll") for (std::uint32_t q = 1; q < ndofs; ++q) N = fj == q ? M_surf[q] : N;
        const double a = AddFaceCouplingToEntry(Ai[c], mu, N, -mu * N);
        Ai[c] = fj < num_face_nodes ? a : Ai[c];
      }
      _Pragma("unroll") for (std::uint32_t fj = 0; fj < ndofs; ++fj)
      {
        const double b = fma(upwind_psi[fj], -mu * M_surf[fj], psi[i]);
        psi[i] = fj < num_face_nodes ? b : psi[i];
      }
    }
    return;
  }
}

/// Compute the sweep matrix from surface integral
template <std::size_t ndofs,
          SweepKind k,
          bool register_resident = (ndofs <= max_register_resident_dofs)>
__CRB_DEVICE_FUNC__ void
ComputeSurfaceIntegral(double* sweep_matrix,
                       double* psi,
                       CellView& cell,
                       DirectionView& direction,
                       const std::uint64_t* cell_edge_data,
                       const unsigned int& angle_group_idx,
                       const unsigned int& group_idx,
                       const Arguments<k>& args)
{
  // loop over each face
  std::uint32_t face_node_counter = 0;
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    // get face view
    FaceView face;
    cell.GetFaceView(face, f);
    // determine if this face is incoming
    NodeIndexType<k> idx(cell_edge_data[face_node_counter]);
    if (not idx.IsUndefined() and not idx.IsOutgoing())
    {
      if constexpr (k == SweepKind::CBC and register_resident and
                    ndofs > max_register_resident_dofs)
        AddIncomingFaceCBC<ndofs, register_resident>(sweep_matrix,
                                                     psi,
                                                     face,
                                                     CosineRoundingSecondProduct(direction, face),
                                                     cell_edge_data + face_node_counter,
                                                     angle_group_idx,
                                                     group_idx,
                                                     args);
      else
      {
        double mu = direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
                    direction.omega[2] * face.normal[2];
        auto add_face_row = [&](std::uint32_t fi, double* Ai, double& psi_i)
        {
          for (std::uint32_t fj = 0; fj < face.num_face_nodes; ++fj)
          {
            std::uint32_t j = face.cell_mapping_data[fj];
            double mu_Nij = -mu * face.M_surf_data[fi * face.num_face_nodes + fj];
            AddAtIndex<ndofs>(Ai, j, mu_Nij);
            double* upwind_psi;
            if constexpr (k == SweepKind::AAH)
              upwind_psi =
                args.flud_data.GetIncomingFluxPointer(cell_edge_data[face_node_counter + fj],
                                                      angle_group_idx,
                                                      group_idx,
                                                      args.boundary,
                                                      args.boundary_offset,
                                                      args.is_surface_source_active);
            else if constexpr (k == SweepKind::CBC)
              upwind_psi =
                args.flud_data.GetIncomingFluxPointer(cell_edge_data[face_node_counter + fj],
                                                      angle_group_idx,
                                                      group_idx,
                                                      args.boundary,
                                                      args.is_surface_source_active);
            psi_i += *upwind_psi * mu_Nij;
          }
        };
        if constexpr (ndofs <= max_register_resident_dofs)
        {
          // Use compile-time row indices, preserving face-node accumulation order in each row.
          _Pragma("unroll") for (std::uint32_t i = 0; i < ndofs; ++i)
          {
            for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
            {
              if (face.cell_mapping_data[fi] == i)
                add_face_row(fi, sweep_matrix + i * ndofs, psi[i]);
            }
          }
        }
        else
        {
          for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
          {
            std::uint32_t i = face.cell_mapping_data[fi];
            add_face_row(fi, sweep_matrix + i * ndofs, psi[i]);
          }
        }
      }
    }
    // update face node counter
    face_node_counter += face.num_face_nodes;
  }
}

/// Call the same device callable as an lvalue with compile-time indices in [begin, end).
template <std::uint32_t begin, std::uint32_t end, class F>
__CRB_DEVICE_FUNC__ void
StaticFor(F&& f) // NOLINT(cppcoreguidelines-missing-std-forward)
{
  if constexpr (begin < end)
  {
    f(std::integral_constant<std::uint32_t, begin>{});
    StaticFor<begin + 1, end>(f);
  }
}

/// Gaussian elimination
template <std::size_t ndofs, bool register_resident = (ndofs <= max_register_resident_dofs)>
__CRB_DEVICE_FUNC__ void
GaussianElimination(double* sweep_matrix, double* psi)
{
  if constexpr (register_resident and ndofs > max_register_resident_dofs)
  {
    // Preserve elimination order with compile-time indices so the system can stay in registers.
    constexpr std::uint32_t n = ndofs;
    double* A = sweep_matrix;
    StaticFor<0, n>(
      [&](auto ic)
      {
        constexpr std::uint32_t i = decltype(ic)::value;
        double inv_diag = 1.0 / A[i * n + i];
        StaticFor<i, n>([&](auto jc) { A[i * n + decltype(jc)::value] *= inv_diag; });
        psi[i] *= inv_diag;
        StaticFor<i + 1, n>(
          [&](auto kc)
          {
            constexpr std::uint32_t r = decltype(kc)::value;
            double factor = -A[r * n + i];
            StaticFor<i, n>(
              [&](auto jc)
              {
                constexpr std::uint32_t j = decltype(jc)::value;
                A[r * n + j] += factor * A[i * n + j];
              });
            psi[r] += factor * psi[i];
          });
      });
    StaticFor<0, n - 1>(
      [&](auto rc)
      {
        constexpr std::uint32_t j = n - 2 - decltype(rc)::value;
        StaticFor<j + 1, n>(
          [&](auto ic)
          {
            constexpr std::uint32_t i = decltype(ic)::value;
            psi[j] -= A[j * n + i] * psi[i];
          });
      });
    return;
  }
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
template <std::size_t ndofs,
          SweepKind k,
          bool register_resident = (ndofs <= max_register_resident_dofs)>
__CRB_DEVICE_FUNC__ void
WritePsiToFludsAndOutflow(const double* psi,
                          CellView& cell,
                          DirectionView& direction,
                          const std::uint64_t* cell_edge_data,
                          const unsigned int& angle_group_idx,
                          const unsigned int& group_idx,
                          const Arguments<k>& args)
{
  // loop over each face
  std::uint32_t face_node_counter = 0;
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    // get face view
    FaceView face;
    cell.GetFaceView(face, f);
    // determine if this face is outgoing
    NodeIndexType<k> idx(cell_edge_data[face_node_counter]);
    if (not idx.IsUndefined() and idx.IsOutgoing())
    {
      // compute outflow for boundary faces
      if (face.outflow != nullptr)
      {
        double mu;
        if constexpr (k == SweepKind::CBC and register_resident and
                      ndofs > max_register_resident_dofs)
          mu = CosineRoundingSecondProduct(direction, face);
        else
          mu = direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
               direction.omega[2] * face.normal[2];
        for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
        {
          const std::uint32_t i = face.cell_mapping_data[fi];
          const double outflow = direction.weight * mu * face.IntS_shapeI_data[fi] *
                                 ValueAtIndex<ndofs, register_resident>(psi, i);
          crb::atomic_add(face.outflow + args.groupset_start + group_idx, outflow);
        }
      }
      // skip for non reflecting boundary
      bool write_downwind_psi = true;
      if constexpr (k == SweepKind::AAH)
      {
        if (idx.IsBoundary() && !idx.IsReflecting())
        {
          write_downwind_psi = false;
        }
      }
      else if constexpr (k == SweepKind::CBC)
        write_downwind_psi = args.flud_data.ShouldStoreOutgoingFlux(idx);
      // write downwind psi
      if (write_downwind_psi)
      {
        for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
        {
          std::uint32_t i = face.cell_mapping_data[fi];
          NodeIndexType<k> node_index(cell_edge_data[face_node_counter + fi]);
          // put copy psi to FLUDS
          double* downwind_psi;
          if constexpr (k == SweepKind::AAH)
            downwind_psi = args.flud_data.GetOutgoingFluxPointer(
              node_index, angle_group_idx, args.boundary, args.boundary_offset);
          else if constexpr (k == SweepKind::CBC)
            downwind_psi =
              args.flud_data.GetOutgoingFluxPointer(node_index, angle_group_idx, args.boundary);
          *downwind_psi = ValueAtIndex<ndofs, register_resident>(psi, i);
        }
      }
    }
    face_node_counter += face.num_face_nodes;
  }
}

/// Compute the scalar flux.
template <std::size_t ndofs, SweepKind k>
__CRB_DEVICE_FUNC__ void
ComputePhi(const double* psi,
           CellView& cell,
           DirectionView& direction,
           const unsigned int& group_idx,
           const std::uint32_t& num_moments,
           const Arguments<k>& args)
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
__CRB_DEVICE_FUNC__ void
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

/// Use compile-time indexing for register-resident local systems, including fixed-DOF CBC.

template <std::size_t ndofs,
          SweepKind k,
          bool register_resident = (ndofs <= max_register_resident_dofs)>
__CRB_DEVICE_FUNC__ void
Sweep(const Arguments<k>& args,
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
  ComputeGMS<ndofs, k>(
    buffer.A(), buffer.b(), buffer.s(), cell, direction, group_idx, num_moments, args);
  ComputeSurfaceIntegral<ndofs, k, register_resident>(
    buffer.A(), buffer.b(), cell, direction, cell_edge_data, angle_group_idx, group_idx, args);
  // solve for the angular flux
  GaussianElimination<ndofs, register_resident>(buffer.A(), buffer.b());
  // save the result
  WritePsiToFludsAndOutflow<ndofs, k, register_resident>(
    buffer.b(), cell, direction, cell_edge_data, angle_group_idx, group_idx, args);
  ComputePhi<ndofs, k>(buffer.b(), cell, direction, group_idx, num_moments, args);
  if (saved_psi != nullptr)
  {
    SaveAngularFlux<ndofs>(
      buffer.b(), saved_psi, cell, angle_group_idx, args.flud_data.stride_size);
  }
}

// NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast, modernize-use-auto)

} // namespace opensn::gpu_kernel
