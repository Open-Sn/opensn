// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/gpu/device_binding.h"
#include "caribou/main.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace crb = caribou;

namespace opensn
{

namespace
{

// MIG instances of one GPU share its PCI bus ID but have their own UUIDs
std::vector<std::string>
GetGPUHardwareIDs(std::uint32_t num_gpus)
{
  std::vector<std::string> ids;
  ids.reserve(num_gpus);
  for (std::uint32_t device = 0; device < num_gpus; ++device)
    ids.push_back(crb::get_pci_bus_id(device) + '/' + crb::get_uuid(device));
  return ids;
}

} // namespace

std::vector<std::string>
GetVisibleGPUIDs()
{
  const std::uint32_t num_gpus = crb::get_num_gpus();
#if defined(__NVCC__) || defined(__HIPCC__)
  return GetGPUHardwareIDs(num_gpus);
#else
  // The PCI bus ID and UUID queries of SYCL are extensions that a device may not support. Device
  // indices then stand in for them. They identify the GPUs correctly when every rank sees all the
  // GPUs of its node, or when each rank sees only its own.
  try
  {
    return GetGPUHardwareIDs(num_gpus);
  }
  catch (const sycl::exception& error)
  {
    if (error.code() != sycl::errc::feature_not_supported)
      throw;
  }
  std::vector<std::string> ids;
  ids.reserve(num_gpus);
  for (std::uint32_t device = 0; device < num_gpus; ++device)
    ids.push_back("index/" + std::to_string(device));
  return ids;
#endif
}

void
SetCurrentGPU(int device)
{
  crb::set_current_device(static_cast<std::uint32_t>(device));
}

bool
GPUBackendIsCUDAOrHIP()
{
#if defined(__NVCC__) || defined(__HIPCC__)
  return true;
#else
  return false;
#endif
}

} // namespace opensn
