// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/debug_tool/print_quadrature.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/quadrature_view.h"
#include <cstddef>

namespace opensn
{

namespace device_dbg
{

__global__ void
PrintQuadrature(char* quad_data)
{
  QuadratureView quad(quad_data);
  std::printf("[DB] Quadrature(num_angles: %" PRIu32 ", num_moments: %" PRIu32 ")\n",
              quad.num_angles,
              quad.num_moments);
  for (std::uint32_t angle_idx = 0; angle_idx < quad.num_angles; ++angle_idx)
  {
    DirectionView direction;
    quad.GetDirectionView(direction, angle_idx);
    std::printf("[DB]   Direction %" PRIu32 ":\n", angle_idx);
    std::printf(
      "[DB]     Omega: %f %f %f\n", direction.omega[0], direction.omega[1], direction.omega[2]);
    std::printf("[DB]     weight: %f\n", direction.weight);
    std::printf("[DB]     m2d: [");
    for (std::uint32_t m = 0; m < quad.num_moments; ++m)
    {
      std::printf("%s%f", ((m != 0) ? " " : ""), direction.m2d[m]);
    }
    std::printf("]\n");
    std::printf("[DB]     d2m: [");
    for (std::uint32_t m = 0; m < quad.num_moments; ++m)
    {
      std::printf("%s%f", ((m != 0) ? " " : ""), direction.d2m[m]);
    }
    std::printf("]\n");
  }
}

} // namespace device_dbg

void
device_dbg::PrintQuadratureOnDevice(QuadratureCarrier& quad)
{
  device_dbg::PrintQuadrature<<<1, 1>>>(quad.GetDevicePtr());
}

} // namespace opensn
