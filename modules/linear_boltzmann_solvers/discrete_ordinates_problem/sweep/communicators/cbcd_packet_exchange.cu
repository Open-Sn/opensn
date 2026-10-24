// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_packet_exchange.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"

namespace opensn::cbcd_packet
{

namespace
{

/// Copy each descriptor's values, with one block per descriptor.
__CRB_GLOBAL_FUNC__ void
CBCDDeviceCopyKernel(const DeviceCopyDescriptor* copies, std::size_t num_copies)
{
#if defined(__NVCC__) || defined(__HIPCC__)
  const auto thread = threadIdx.x;
  const auto threads = blockDim.x;
  const auto first_face = blockIdx.x;
  const auto face_stride = gridDim.x;
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  auto work_index = ::sycl::ext::oneapi::this_work_item::get_nd_item<3>();
  const auto thread = work_index.get_local_id(2);
  const auto threads = work_index.get_local_range(2);
  const auto first_face = work_index.get_group(2);
  const auto face_stride = work_index.get_group_range(2);
#endif
  for (std::size_t face = first_face; face < num_copies; face += face_stride)
  {
    const auto copy = copies[face];
    for (std::size_t i = thread; i < copy.num_values; i += threads)
      copy.destination[i] = copy.source[i];
  }
}

/// Write each view's node offsets, with one thread per view.
__CRB_GLOBAL_FUNC__ void
CBCDPacketViewKernel(const PacketView* views, std::size_t num_views)
{
#if defined(__NVCC__) || defined(__HIPCC__)
  const std::size_t first_face = threadIdx.x + blockIdx.x * blockDim.x;
  const std::size_t face_stride = blockDim.x * gridDim.x;
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  auto work_index = ::sycl::ext::oneapi::this_work_item::get_nd_item<3>();
  const auto first_face = work_index.get_global_id(2);
  const auto face_stride = work_index.get_global_range(2);
#endif
  for (std::size_t face = first_face; face < num_views; face += face_stride)
  {
    const auto view = views[face];
    for (std::size_t node = 0; node < view.num_nodes; ++node)
      view.node_offsets[node] = view.packet_offset + node * view.stride;
  }
}

} // namespace

void
DeviceCopyBatch::Submit(crb::Stream& stream)
{
  if (num_copies_ == 0)
    return;
  const auto num_copies = num_copies_;
  auto* copies = device_descriptors_.get();
  crb::impl::copy_h2d(copies, descriptors_.data(), num_copies, stream);
  crb::Dim3 grid_size(static_cast<unsigned int>(std::min<std::size_t>(num_copies, 65535)));
  crb::Dim3 block_size(gpu_kernel::threshold);
#if defined(__NVCC__) || defined(__HIPCC__)
  CBCDDeviceCopyKernel<<<grid_size, block_size, 0, stream>>>(copies, num_copies);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  stream.parallel_for(sycl::nd_range<3>(grid_size * block_size, block_size),
                      [=](sycl::nd_item<3> work_index)
                      { CBCDDeviceCopyKernel(copies, num_copies); });
#endif
}

void
PacketViewBatch::Submit(crb::Stream& stream)
{
  if (num_views_ == 0)
    return;
  const auto num_views = num_views_;
  auto* views = device_views_.get();
  crb::impl::copy_h2d(views, views_.data(), num_views, stream);
  const auto num_blocks = (num_views + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
  crb::Dim3 grid_size(static_cast<unsigned int>(std::min<std::size_t>(num_blocks, 65535)));
  crb::Dim3 block_size(gpu_kernel::threshold);
#if defined(__NVCC__) || defined(__HIPCC__)
  CBCDPacketViewKernel<<<grid_size, block_size, 0, stream>>>(views, num_views);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  stream.parallel_for(sycl::nd_range<3>(grid_size * block_size, block_size),
                      [=](sycl::nd_item<3> work_index) { CBCDPacketViewKernel(views, num_views); });
#endif
}

} // namespace opensn::cbcd_packet
