#include "modules/linear_boltzmann_solvers/a_lbs_solver/source_functions/adjoint_src_function.h"

#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/lbs_adj_solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace lbs
{

AdjointSourceFunction::AdjointSourceFunction(const LBSSolver& lbs_solver)
  : SourceFunction(lbs_solver)
{
}

void
AdjointSourceFunction::AddVolumetricQOISources(LBSGroupset& groupset,
                                               std::vector<double>& destination_q,
                                               const std::vector<double>&,
                                               int source_flags)
{
  const std::string fname = "AdjointSourceFunction::AddVolumetricQOISources";

  typedef const DiscreteOrdinatesAdjointSolver CAdjointSolver;
  const auto adjoint_solver_ptr = dynamic_cast<CAdjointSolver*>(&lbs_solver_);
  if (not adjoint_solver_ptr)
    throw std::logic_error(fname + ": Failed to cast lbs_solver_ to adjoint.");

  const auto& adjoint_solver = *adjoint_solver_ptr;

#ifdef OPENSN_WITH_LUA
  const auto& response_functions = adjoint_solver.GetResponseFunctions();
  const auto& basic_options = adjoint_solver.GetBasicOptions();
  const auto& cell_transport_views = adjoint_solver.GetCellTransportViews();
  const auto& grid = adjoint_solver.Grid();
  const auto num_groups = adjoint_solver.NumGroups();

  const bool apply_fixed_src = (source_flags & APPLY_FIXED_SOURCES);

  const auto gs_i = static_cast<size_t>(groupset.groups_.front().id_);
  const auto gs_f = static_cast<size_t>(groupset.groups_.back().id_);

  if (apply_fixed_src)
    for (const auto& qoi_data : response_functions)
    {
      const auto& qoi_designation = qoi_data.first;
      const auto& qoi_cell_subscription = qoi_data.second;

      if (qoi_designation.name == basic_options("REFERENCE_RF").StringValue())
      {
        for (size_t local_id : qoi_cell_subscription)
        {
          const auto& full_cell_view = cell_transport_views[local_id];
          const auto& cell = grid.local_cells[local_id];
          const auto& response = qoi_designation.GetMGResponse(cell, num_groups);
          const int num_nodes = full_cell_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            size_t uk_map = full_cell_view.MapDOF(i, 0, 0); // unknown map
            for (size_t g = gs_i; g <= gs_f; ++g)
              destination_q[uk_map + g] += response[g];
          } // for node
        }   // for local cell-id of qoi
      }     // if ref-qoi
    }       // for qoi
#endif
}

} // namespace lbs
