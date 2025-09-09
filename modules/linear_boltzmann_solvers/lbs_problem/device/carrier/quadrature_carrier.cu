// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"

namespace opensn
{

QuadratureCarrier::QuadratureCarrier(const LBSGroupset& groupset)
{
  std::uint64_t size = ComputeSize(groupset);
  host_memory_.reserve(size);
  host_memory_.resize(size);
  Assemble(groupset);
  device_memory_ = crb::DeviceMemory<char>(size);
  crb::copy(device_memory_, host_memory_, size);
}

std::uint64_t
QuadratureCarrier::ComputeSize(const LBSGroupset& groupset)
{
  std::uint64_t alloc_size = 0;
  // number of angles and number of moments
  alloc_size += 2 * sizeof(std::uint32_t);
  // get number of angles and number of moments
  const AngularQuadrature& quadrature = *(groupset.quadrature);
  std::size_t num_angles = quadrature.omegas.size();
  std::size_t num_moments = quadrature.GetMomentToDiscreteOperator().size();
  // size of each directions
  alloc_size += num_angles * (4 * sizeof(double) + 2 * num_moments * sizeof(double));
  return alloc_size;
}

void
QuadratureCarrier::Assemble(const LBSGroupset& groupset)
{
  // get information
  const AngularQuadrature& quadrature = *(groupset.quadrature);
  // number of angles and number of moments
  char* data = reinterpret_cast<char*>(host_memory_.data());
  std::uint32_t* num_angles_and_moments_data = reinterpret_cast<std::uint32_t*>(data);
  std::size_t num_angles = quadrature.omegas.size();
  *(num_angles_and_moments_data++) = num_angles;
  const std::vector<std::vector<double>>& m2d = quadrature.GetMomentToDiscreteOperator();
  const std::vector<std::vector<double>>& d2m = quadrature.GetDiscreteToMomentOperator();
  std::size_t num_moments = quadrature.GetNumMoments();
  *(num_angles_and_moments_data++) = num_moments;
  data = reinterpret_cast<char*>(num_angles_and_moments_data);
  // direction data
  for (std::size_t direction_num = 0; direction_num < num_angles; ++direction_num)
  {
    // omega
    const Vector3& omega = quadrature.omegas[direction_num];
    double* omega_data = reinterpret_cast<double*>(data);
    *(omega_data++) = omega.x;
    *(omega_data++) = omega.y;
    *(omega_data++) = omega.z;
    // weight
    double* weight_data = omega_data;
    *(weight_data++) = quadrature.weights[direction_num];
    // M2D data
    double* m2d_data = weight_data;
    for (std::size_t m = 0; m < num_moments; ++m)
    {
      *(m2d_data++) = m2d[direction_num][m];
    }
    // D2M data
    double* d2m_data = m2d_data;
    for (std::size_t m = 0; m < num_moments; ++m)
    {
      *(d2m_data++) = d2m[direction_num][m];
    }
    data = reinterpret_cast<char*>(d2m_data);
  }
}

} // namespace opensn
