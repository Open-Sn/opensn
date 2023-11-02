#include "CBC_SPDS.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/SweepUtilities/sweep_namespace.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "graphs/chi_directed_graph.h"
#include "utils/chi_timer.h"

namespace lbs
{

CBC_SPDS::CBC_SPDS(const chi_mesh::Vector3& omega,
                   const chi_mesh::MeshContinuum& grid,
                   bool cycle_allowance_flag,
                   bool verbose)
  : SPDS(omega, grid, verbose)
{
  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Building sweep ordering for Omega = "
                          << omega.PrintS();

  size_t num_loc_cells = grid.local_cells.size();

  //============================================= Populate Cell Relationships
  Chi::log.Log0Verbose1() << "Populating cell relationships";
  std::vector<std::set<std::pair<int, double>>> cell_successors(num_loc_cells);
  std::set<int> location_successors;
  std::set<int> location_dependencies;

  PopulateCellRelationships(
    omega, location_dependencies, location_successors, cell_successors);

  location_successors_.reserve(location_successors.size());
  location_dependencies_.reserve(location_dependencies.size());

  for (auto v : location_successors)
    location_successors_.push_back(v);

  for (auto v : location_dependencies)
    location_dependencies_.push_back(v);

  //============================================= Build graph
  chi::DirectedGraph local_DG;

  // Add vertex for each local cell
  for (int c = 0; c < num_loc_cells; ++c)
    local_DG.AddVertex();

  // Create graph edges
  for (int c = 0; c < num_loc_cells; c++)
    for (auto& successor : cell_successors[c])
      local_DG.AddEdge(c, successor.first, successor.second);

  //============================================= Remove local cycles if allowed
  if (verbose_) PrintedGhostedGraph();

  if (cycle_allowance_flag)
  {
    Chi::log.Log0Verbose1()
      << Chi::program_timer.GetTimeString() << " Removing inter-cell cycles.";

    auto edges_to_remove = local_DG.RemoveCyclicDependencies();

    for (auto& edge_to_remove : edges_to_remove)
    {
      local_cyclic_dependencies_.emplace_back(edge_to_remove.first,
                                              edge_to_remove.second);
    }
  }

  //============================================= Generate topological sorting
  Chi::log.Log0Verbose1()
    << Chi::program_timer.GetTimeString()
    << " Generating topological sorting for local sweep ordering";
  auto so_temp = local_DG.GenerateTopologicalSort();
  spls_.item_id.clear();
  for (auto v : so_temp)
    spls_.item_id.emplace_back(v);

  if (spls_.item_id.empty())
  {
    Chi::log.LogAllError()
      << "Topological sorting for local sweep-ordering failed. "
      << "Cyclic dependencies detected. Cycles need to be allowed"
      << " by calling application.";
    Chi::Exit(EXIT_FAILURE);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Task
  //                                                        Dependency Graphs
  // All locations will gather other locations' dependencies
  // so that each location has the ability to build
  // the global task graph.

  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Communicating sweep dependencies.";

  // auto& global_dependencies = sweep_order->global_dependencies;
  std::vector<std::vector<int>> global_dependencies;
  global_dependencies.resize(Chi::mpi.process_count);

  chi_mesh::sweep_management::CommunicateLocationDependencies(
    location_dependencies_, global_dependencies);

  ////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Build task
  ////                                                        dependency graph
  // BuildTaskDependencyGraph(global_dependencies, cycle_allowance_flag);

  constexpr auto INCOMING =
    chi_mesh::sweep_management::FaceOrientation::INCOMING;
  constexpr auto OUTGOING =
    chi_mesh::sweep_management::FaceOrientation::OUTGOING;

  // For each local cell create a task
  for (const auto& cell : grid_.local_cells)
  {
    const size_t num_faces = cell.faces_.size();
    unsigned int num_dependencies = 0;
    std::vector<uint64_t> succesors;

    for (size_t f = 0; f < num_faces; ++f)
      if (cell_face_orientations_[cell.local_id_][f] == INCOMING)
      {
        if (cell.faces_[f].has_neighbor_)
          ++num_dependencies;

      }
      else if (cell_face_orientations_[cell.local_id_][f] == OUTGOING)
      {
        const auto& face = cell.faces_[f];
        if (face.has_neighbor_ and grid.IsCellLocal(face.neighbor_id_))
          succesors.push_back(grid.cells[face.neighbor_id_].local_id_);
      }

    task_list_.push_back({num_dependencies,
                          succesors,
                          /*reference_id_=*/cell.local_id_,
                          /*cell_ptr_=*/&cell,
                          /*completed_=*/false});
  } // for cell in SPLS

  Chi::mpi.Barrier();

  Chi::log.Log0Verbose1() << Chi::program_timer.GetTimeString()
                          << " Done computing sweep ordering.\n\n";
}

const std::vector<chi_mesh::sweep_management::Task>& CBC_SPDS::TaskList() const
{
  return task_list_;
}

} // namespace lbs