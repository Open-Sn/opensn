// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_sweep_launch.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_grid_stride.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"
#include "caribou/main.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <type_traits>

// Read launch metadata without a private copy per GPU thread.
#ifdef __NVCC__
#define CBCDGridConstant __grid_constant__
#else
#define CBCDGridConstant
#endif

#ifdef __HIPCC__
#define CBCDLaunchBounds(fixed)                                                                    \
  __attribute__((                                                                                  \
    amdgpu_flat_work_group_size(opensn::gpu_kernel::CBCD_HIP_MIN_BLOCK_THREADS,                    \
                                (fixed) > opensn::gpu_kernel::max_register_resident_dofs           \
                                  ? opensn::gpu_kernel::threshold                                  \
                                  : opensn::gpu_kernel::CBCD_HIP_MAX_BLOCK_THREADS)))
#else
#define CBCDLaunchBounds(fixed)
#endif

namespace crb = caribou;

namespace opensn::gpu_kernel
{

/// Sentinel selecting the kernel that dispatches every supported cell DOF count.
inline constexpr std::uint32_t CBCD_ALL_CELL_DOFS = 0;

/// HIP workgroup bound for the default register policy. Larger fixed-DOF systems use the
/// dispatched block size.
inline constexpr unsigned int CBCD_HIP_MAX_BLOCK_THREADS = 1024;
inline constexpr unsigned int CBCD_DYNAMIC_SHARED_BYTES = 0;

inline constexpr unsigned int CBCD_HIP_MIN_BLOCK_THREADS = 1;

/// Explicit kernel argument layout for size and alignment checks.
/// CUDA and HIP pass metadata by value. SYCL passes a device pointer.
template <typename Launch>
struct CBCD_ExplicitArguments
{
#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  const Launch* launch;
#else
  Launch launch;
#endif
  const Arguments<SweepKind::CBC>* arguments;
  double* const* saved_psi;
  std::uint64_t num_blocks;
};

template <typename Launch>
constexpr bool
IsCBCDExplicitArgumentLayoutValid()
{
  using Explicit = CBCD_ExplicitArguments<Launch>;
  return std::is_standard_layout_v<Explicit> and std::is_trivially_copyable_v<Explicit> and
         offsetof(Explicit, launch) == 0 and
         offsetof(Explicit, arguments) == sizeof(Explicit::launch) and
         sizeof(Explicit) - sizeof(Explicit::launch) <= CBCD_OTHER_PARAMETER_ALLOWANCE;
}

static_assert(IsCBCDExplicitArgumentLayoutValid<CBCD_LaunchClass<0>>() and
                IsCBCDExplicitArgumentLayoutValid<CBCD_LaunchClass<1>>() and
                IsCBCDExplicitArgumentLayoutValid<CBCD_LaunchClass<2>>(),
              "Device CBC's explicit kernel arguments exceed their selected allowance.");

/// Borrow launch metadata for this kernel invocation. It must be the first kernel argument
/// so AMD GPUs can read it at offset zero in the kernarg segment.
template <typename Launch>
__CRB_DEVICE_FUNC__ const Launch*
GetKernelLaunch(const Launch* launch)
{
  static_assert(std::is_trivially_copyable_v<Launch> and std::is_standard_layout_v<Launch>);
#if defined(__HIP_DEVICE_COMPILE__) && defined(__AMDGCN__)
  (void)launch;
  return reinterpret_cast<const Launch*>(__builtin_amdgcn_kernarg_segment_ptr());
#else
  return launch;
#endif
}

/// Sweep one logical block of a fused launch using its read-only metadata.
template <typename Launch, std::uint32_t fixed_ndofs>
__CRB_DEVICE_FUNC__ void
CBCDSweepBlock(const Launch& launch,
               const Arguments<SweepKind::CBC>* __restrict__ arguments,
               double* const* __restrict__ saved_psi,
               const std::uint64_t block_idx,
               const unsigned int thread_idx)
{
  std::uint32_t begin = 0;
  std::uint32_t end = launch.num_slots;
  while (begin < end)
  {
    const auto mid = begin + (end - begin) / 2;
    if (block_idx < launch.slots[mid].block_end)
      end = mid;
    else
      begin = mid + 1;
  }
  const auto& slot = launch.slots[begin];
  const std::uint32_t block_size_x = slot.block_size_x;
  const std::uint32_t num_stride_blocks = slot.num_stride_blocks;
  // RoundUp can produce a block width that does not divide threshold.
  if (thread_idx >= block_size_x * (threshold / block_size_x))
    return;
  const auto local_block_idx = block_idx - (begin == 0 ? 0 : launch.slots[begin - 1].block_end);
  auto args = arguments[slot.angle_set_id];
  const auto cell_idx =
    thread_idx / block_size_x + (threshold / block_size_x) * (local_block_idx / num_stride_blocks);
  const auto angle_group_idx = static_cast<unsigned int>(
    thread_idx % block_size_x + block_size_x * (local_block_idx % num_stride_blocks));
  if (cell_idx >= slot.num_cells or angle_group_idx >= args.flud_data.stride_size)
    return;

  const auto angle_idx = angle_group_idx / args.groupset_size;
  auto group_idx = angle_group_idx - angle_idx * args.groupset_size;
  const auto cell_data = launch.entries[slot.first_entry + cell_idx];
  const auto cell_local_idx = static_cast<std::uint32_t>(cell_data);
  const auto* faces = &launch.entries[slot.first_entry + slot.num_cells];
  args.flud_data.incoming_payload_offsets =
    faces + ((cell_data >> CBCD_CELL_ID_BITS) & CBCD_FACE_OFFSET_MASK);
  args.flud_data.outgoing_payload_offsets =
    faces + slot.num_incoming_faces + (cell_data >> CBCD_OUTGOING_OFFSET_SHIFT);
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  if constexpr (fixed_ndofs == CBCD_ALL_CELL_DOFS)
  {
    if (cell.num_nodes == 0)
      return;
  }
  else if (cell.num_nodes != fixed_ndofs)
    return;
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  QuadratureView quadrature(args.quad_data);
  auto num_moments = quadrature.num_moments;
  DirectionView direction;
  quadrature.GetDirectionView(direction, args.directions[angle_idx]);
  auto* slot_saved_psi = saved_psi[slot.angle_set_id];
  if constexpr (fixed_ndofs == CBCD_ALL_CELL_DOFS)
  {
    auto stride_idx = angle_group_idx;
    SweepDispatch<SweepKind::CBC>(cell.num_nodes,
                                  args,
                                  cell,
                                  direction,
                                  cell_edge_data,
                                  stride_idx,
                                  group_idx,
                                  num_moments,
                                  slot_saved_psi);
  }
  else
    Sweep<fixed_ndofs, SweepKind::CBC, true>(args,
                                             cell,
                                             direction,
                                             cell_edge_data,
                                             angle_group_idx,
                                             group_idx,
                                             num_moments,
                                             slot_saved_psi);
}

/// Sweep ready cells. Each logical block belongs to one slot and covers cells and a tile of
/// angle/group pairs. SYCL strides over logical blocks. CUDA and HIP handle one per grid block,
/// letting compilers load unmodified sweep arguments before any stores.
///
/// The argument table and the saved-flux pointer table are separate read-only allocations.
/// fixed_ndofs selects a register-resident system for that cell size. Zero dispatches per cell
/// with the default register policy. Only dispatch a fixed count on a matching local mesh.
#if defined(__NVCC__) || defined(__HIPCC__)
template <typename Launch, std::uint32_t fixed_ndofs>
__CRB_GLOBAL_FUNC__ void
CBCDFusedSweepKernel(const CBCDGridConstant Launch launch_parameter,
                     const Arguments<SweepKind::CBC>* __restrict__ arguments,
                     double* const* __restrict__ saved_psi,
                     std::uint64_t num_blocks) CBCDLaunchBounds(fixed_ndofs)
{
  const Launch& launch = *GetKernelLaunch(&launch_parameter);
#if defined(__NVCC__) || defined(__HIPCC__)
  static_cast<void>(num_blocks);
  CBCDSweepBlock<Launch, fixed_ndofs>(launch, arguments, saved_psi, blockIdx.x, threadIdx.x);
#endif
}

#endif

#undef CBCDGridConstant
#undef CBCDLaunchBounds

/// Enqueue on the owning worker's stream. CUDA and HIP copy launch metadata on submission.
/// On SYCL the lane retains device metadata through queue completion. Keep the device argument
/// and pointer tables and the buffers they reference alive through completion on every backend.
template <typename Launch, std::uint32_t fixed_ndofs>
void
LaunchCBCDFusedSweep(const Launch& launch,
                     const Arguments<SweepKind::CBC>* arguments,
                     double* const* saved_psi,
                     CBCD_LaunchStream& stream)
{
  assert(launch.num_slots != 0);
  const std::uint64_t num_blocks = launch.slots[launch.num_slots - 1].block_end;
#if defined(__NVCC__) || defined(__HIPCC__)
  crb::Dim3 grid_size(static_cast<unsigned int>(num_blocks));
#else
  crb::Dim3 grid_size(
    static_cast<unsigned int>(std::min<std::uint64_t>(num_blocks, stream.grid_blocks)));
#endif
  crb::Dim3 block_size(threshold);
#if defined(__NVCC__) || defined(__HIPCC__)
  CBCDFusedSweepKernel<Launch, fixed_ndofs>
    <<<grid_size, block_size, CBCD_DYNAMIC_SHARED_BYTES, stream>>>(
      launch, arguments, saved_psi, num_blocks);
  crb::impl::check_error(::GPU_API(GetLastError)());
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  static_assert(std::is_same_v<Launch, CBCD_Launch>);
  const auto* metadata = stream.metadata;
  stream.queue.parallel_for(sycl::nd_range<3>(grid_size * block_size, block_size),
                            [=](sycl::nd_item<3> work_index)
                            {
                              CBCDGridStride(
                                work_index,
                                num_blocks,
                                [=](std::uint64_t block_idx, unsigned int thread_idx)
                                {
                                  CBCDSweepBlock<Launch, fixed_ndofs>(
                                    *metadata, arguments, saved_psi, block_idx, thread_idx);
                                });
                            });
#endif
}

/// Launch using a size class that fits the metadata. Keep each all-counts kernel in its own
/// translation unit so sweep dispatch can inline its arguments. Compile fixed-DOF kernels
/// separately to isolate register allocation from the all-counts kernels.
template <std::size_t launch_class, std::uint32_t fixed_ndofs>
void
LaunchCBCDFusedSweepAs(const CBCD_Launch& launch,
                       const Arguments<SweepKind::CBC>* arguments,
                       double* const* saved_psi,
                       CBCD_LaunchStream& stream)
{
  if constexpr (launch_class + 1 == CBCD_NUM_LAUNCH_CLASSES)
    LaunchCBCDFusedSweep<CBCD_Launch, fixed_ndofs>(launch, arguments, saved_psi, stream);
  else if constexpr (IsCBCDLaunchClassUsed(launch_class))
  {
    CBCD_LaunchClass<launch_class> sized;
    sized.CopyUsed(launch);
    LaunchCBCDFusedSweep<CBCD_LaunchClass<launch_class>, fixed_ndofs>(
      sized, arguments, saved_psi, stream);
  }
}

extern template void LaunchCBCDFusedSweepAs<0, CBCD_ALL_CELL_DOFS>(const CBCD_Launch&,
                                                                   const Arguments<SweepKind::CBC>*,
                                                                   double* const*,
                                                                   CBCD_LaunchStream&);
extern template void LaunchCBCDFusedSweepAs<1, CBCD_ALL_CELL_DOFS>(const CBCD_Launch&,
                                                                   const Arguments<SweepKind::CBC>*,
                                                                   double* const*,
                                                                   CBCD_LaunchStream&);

} // namespace opensn::gpu_kernel
