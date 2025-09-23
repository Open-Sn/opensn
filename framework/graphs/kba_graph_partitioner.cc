// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/kba_graph_partitioner.h"
#include "framework/utils/utils.h"
#include "framework/mesh/mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <format>
#include <cmath>

namespace opensn
{

OpenSnRegisterObjectInNamespace(mesh, KBAGraphPartitioner);

InputParameters
KBAGraphPartitioner::GetInputParameters()
{
  InputParameters params = GraphPartitioner::GetInputParameters();

  params.SetGeneralDescription("Koch, Baker and Alcouffe based partitioning. "
                               "This is an overlayed ortho-grid based partitioner");

  params.AddOptionalParameter("nx", 1, "The number of partitions in x");
  params.AddOptionalParameter("ny", 1, "The number of partitions in y");
  params.AddOptionalParameter("nz", 1, "The number of partitions in z");

  params.AddOptionalParameter(
    "xcuts", std::vector<double>{}, "Location of the internal x-cuts. Require nx-1 entries");
  params.AddOptionalParameter(
    "ycuts", std::vector<double>{}, "Location of the internal y-cuts. Require ny-1 entries");
  params.AddOptionalParameter(
    "zcuts", std::vector<double>{}, "Location of the internal z-cuts. Require nz-1 entries");

  return params;
}

std::shared_ptr<KBAGraphPartitioner>
KBAGraphPartitioner::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<KBAGraphPartitioner>("mesh::KBAGraphPartitioner", params);
}

KBAGraphPartitioner::KBAGraphPartitioner(const InputParameters& params)
  : GraphPartitioner(params),
    nx_(params.GetParamValue<size_t>("nx")),
    ny_(params.GetParamValue<size_t>("ny")),
    nz_(params.GetParamValue<size_t>("nz")),
    xcuts_(params.GetParamVectorValue<double>("xcuts")),
    ycuts_(params.GetParamVectorValue<double>("ycuts")),
    zcuts_(params.GetParamVectorValue<double>("zcuts")),
    coordinate_infos_{CoordinateInfo{&xcuts_, nx_, "x"},
                      CoordinateInfo{&ycuts_, ny_, "y"},
                      CoordinateInfo{&zcuts_, nz_, "z"}}
{
  for (const auto& [cuts_ptr, n, name] : coordinate_infos_)
  {
    const auto& cuts = *cuts_ptr;

    // Check number of items
    if (cuts.size() != (n - 1))
    {
      auto err = std::format(
        "The number of cuts supplied for \"{}cuts\" is not equal to n{}-1.", name, name);
      OpenSnInvalidArgument(err);
    }
    if (cuts.empty())
      continue;

    // Check monitonically increasing
    {
      double prev_value = 0.0;
      for (const double cut_value : *cuts_ptr)
      {
        OpenSnInvalidArgumentIf(cut_value != cuts.front() and cut_value <= prev_value,
                                "Parameter \"" + name +
                                  "\" requires monotonically increasing values");
        prev_value = cut_value;
      }
    } // for cut value
  } // for each coordinate
}

std::vector<int>
KBAGraphPartitioner::Partition(const std::vector<std::vector<uint64_t>>& graph,
                               const std::vector<Vector3>& centroids,
                               int number_of_parts)
{
  log.Log0Verbose1() << "Partitioning with KBAGraphPartitioner";

  OpenSnLogicalErrorIf(centroids.size() != graph.size(),
                       "Graph number of entries not equal to centroids' number of entries.");
  const size_t num_cells = graph.size();
  std::vector<size_t> pids(num_cells, 0);
  for (size_t c = 0; c < num_cells; ++c)
  {
    const auto& point = centroids[c];
    // Partitions the point
    std::array<size_t, 3> p_vals = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i)
    {
      const auto& cuts = *coordinate_infos_[i].cuts_;
      const size_t num_cuts = cuts.size();

      size_t p_val;
      bool home_found = false;
      for (size_t j = 0; j < num_cuts; ++j)
        if (cuts[j] > point[i])
        {
          p_val = j;
          home_found = true;
          break;
        }

      p_vals[i] = home_found ? p_val : (coordinate_infos_[i].n_ - 1);
    }

    const auto nx = static_cast<int64_t>(coordinate_infos_[0].n_);
    const auto ny = static_cast<int64_t>(coordinate_infos_[1].n_);

    const auto i = static_cast<int64_t>(p_vals[0]);
    const auto j = static_cast<int64_t>(p_vals[1]);
    const auto k = static_cast<int64_t>(p_vals[2]);

    pids[c] = nx * ny * k + nx * j + i;
  } // for cell c

  if ((nx_ * ny_ * nz_) != number_of_parts)
    log.Log0Warning() << "KBAGraphPartitioner::Partition nx_*ny_*nz_ != number_of_parts";

  const auto pid_subsets = MakeSubSets(nx_ * ny_ * nz_, number_of_parts);

  std::vector<int> real_pids(num_cells, 0);
  for (size_t c = 0; c < num_cells; ++c)
  {
    for (int p = 0; p < number_of_parts; ++p)
    {
      if (pids[c] >= pid_subsets[p].ss_begin and pids[c] <= pid_subsets[p].ss_end)
        real_pids[c] = p;
    }
  }

  log.Log0Verbose1() << "Done partitioning with KBAGraphPartitioner";

  return real_pids;
}

} // namespace opensn
