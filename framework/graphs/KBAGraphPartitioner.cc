#include "opensn/framework/graphs/KBAGraphPartitioner.h"

#include "opensn/framework/ChiObjectFactory.h"
#include "opensn/framework/utils/chi_utils.h"

#include "opensn/framework/mesh/chi_mesh.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

#include <cmath>

namespace chi
{

RegisterChiObject(chi, KBAGraphPartitioner);

InputParameters
KBAGraphPartitioner::GetInputParameters()
{
  InputParameters params = GraphPartitioner::GetInputParameters();

  params.SetGeneralDescription("Koch, Baker and Alcouffe based partitioning. "
                               "This is an overlayed ortho-grid based partitioner");
  params.SetDocGroup("Graphs");

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

    //======================= Check number of items
    if (cuts.size() != (n - 1))
      ChiInvalidArgument("The number of cuts supplied for \"" + name + "cuts\" is not equal to n" +
                         name + "-1.");
    if (cuts.empty()) continue;

    //======================= Check monitonically increasing
    {
      double prev_value = 0.0;
      for (const double cut_value : *cuts_ptr)
      {
        ChiInvalidArgumentIf(cut_value != cuts.front() and cut_value <= prev_value,
                             "Parameter \"" + name + "\" requires monotonically increasing values");
        prev_value = cut_value;
      }
    } // for cut value
  }   // for each coordinate
}

std::vector<int64_t>
KBAGraphPartitioner::Partition(const std::vector<std::vector<uint64_t>>& graph,
                               const std::vector<chi_mesh::Vector3>& centroids,
                               int number_of_parts)
{
  Chi::log.Log0Verbose1() << "Partitioning with KBAGraphPartitioner";

  ChiLogicalErrorIf(centroids.size() != graph.size(),
                    "Graph number of entries not equal to centroids' number of entries.");
  const size_t num_cells = graph.size();
  std::vector<int64_t> pids(num_cells, 0);
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

    const int64_t nx = static_cast<int64_t>(coordinate_infos_[0].n_);
    const int64_t ny = static_cast<int64_t>(coordinate_infos_[1].n_);

    const int64_t i = static_cast<int64_t>(p_vals[0]);
    const int64_t j = static_cast<int64_t>(p_vals[1]);
    const int64_t k = static_cast<int64_t>(p_vals[2]);

    pids[c] = nx * ny * k + nx * j + i;
  } // for cell c

  if ((nx_ * ny_ * nz_) != number_of_parts)
    Chi::log.Log0Warning() << "KBAGraphPartitioner::Partition nx_*ny_*nz_ != number_of_parts";

  const auto pid_subsets = chi::MakeSubSets(nx_ * ny_ * nz_, number_of_parts);

  std::vector<int64_t> real_pids(num_cells, 0);
  for (size_t c = 0; c < num_cells; ++c)
  {
    for (size_t p = 0; p < number_of_parts; ++p)
    {
      if (pids[c] >= pid_subsets[p].ss_begin and pids[c] <= pid_subsets[p].ss_end)
        real_pids[c] = static_cast<int64_t>(p);
    }
  }

  Chi::log.Log0Verbose1() << "Done partitioning with KBAGraphPartitioner";

  return real_pids;
}

} // namespace chi
