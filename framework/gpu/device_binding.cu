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

#if defined(__NVCC__) || defined(__HIPCC__)

std::vector<std::string>
GetVisibleGPUBusIDs()
{
  const std::uint32_t num_gpus = crb::get_num_gpus();
  std::vector<std::string> bus_ids;
  bus_ids.reserve(num_gpus);
  for (std::uint32_t device = 0; device < num_gpus; ++device)
    bus_ids.push_back(crb::get_pci_bus_id(device));
  return bus_ids;
}

void
SetCurrentGPU(int device)
{
  crb::set_device(static_cast<std::uint32_t>(device));
}

#else

// caribou's SYCL backend selects its device when the library loads, so BindDevice binds none.
std::vector<std::string>
GetVisibleGPUBusIDs()
{
  return {};
}

void
SetCurrentGPU(int /*device*/)
{
}

#endif

} // namespace opensn
