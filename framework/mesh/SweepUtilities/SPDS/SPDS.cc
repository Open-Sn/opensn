#include "framework/mesh/SweepUtilities/SPDS/SPDS.h"

#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/mpi/chi_mpi.h"
#include "framework/console/chi_console.h"
#include "framework/utils/chi_timer.h"

#include <algorithm>

/** Given a location J index, maps to a predecessor location.*/
int
chi_mesh::sweep_management::SPDS::MapLocJToPrelocI(int locJ) const
{
  for (int i = 0; i < location_dependencies_.size(); i++)
  {
    if (location_dependencies_[i] == locJ) { return i; }
  }

  for (int i = 0; i < delayed_location_dependencies_.size(); i++)
  {
    if (delayed_location_dependencies_[i] == locJ) { return -(i + 1); }
  }

  Chi::log.LogAllError() << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
  Chi::Exit(EXIT_FAILURE);
  return 0;
}

/** Given a location J index, maps to a dependent location.*/
int
chi_mesh::sweep_management::SPDS::MapLocJToDeplocI(int locJ) const
{
  for (int i = 0; i < location_successors_.size(); i++)
  {
    if (location_successors_[i] == locJ) { return i; }
  }

  Chi::log.LogAllError() << "SPDS Invalid mapping encountered in MapLocJToDeplocI.";
  Chi::Exit(EXIT_FAILURE);
  return 0;
}

/**Populates cell relationships*/
void
chi_mesh::sweep_management::SPDS::PopulateCellRelationships(
  const chi_mesh::Vector3& omega,
  std::set<int>& location_dependencies,
  std::set<int>& location_successors,
  std::vector<std::set<std::pair<int, double>>>& cell_successors)
{
  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations_.assign(grid_.local_cells.size(), {});
  for (auto& cell : grid_.local_cells)
    cell_face_orientations_[cell.local_id_].assign(cell.faces_.size(), FOPARALLEL);

  for (auto& cell : grid_.local_cells)
  {
    size_t f = 0;
    for (auto& face : cell.faces_)
    {
      //======================================= Determine if the face
      //                                        is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal_);

      bool owns_face = true;
      if (face.has_neighbor_ and cell.global_id_ > face.neighbor_id_) owns_face = false;

      if (owns_face)
      {
        // clang-format off
        if (mu > tolerance) orientation = FOOUTGOING;
        else if (mu < tolerance) orientation = FOINCOMING;

        cell_face_orientations_[cell.local_id_][f] = orientation;

        if (face.has_neighbor_ and grid_.IsCellLocal(face.neighbor_id_))
        {
          const auto& adj_cell = grid_.cells[face.neighbor_id_];
          const auto ass_face = face.GetNeighborAssociatedFace(grid_);
          auto& adj_face_ori =
            cell_face_orientations_[adj_cell.local_id_][ass_face];

          switch (orientation)
          {
            case FOPARALLEL: adj_face_ori = FOPARALLEL; break;
            case FOINCOMING: adj_face_ori = FOOUTGOING; break;
            case FOOUTGOING: adj_face_ori = FOINCOMING; break;
          }
        }
        // clang-format on
      } // if face owned
      else if (face.has_neighbor_ and not grid_.IsCellLocal(face.neighbor_id_))
      {
        const auto& adj_cell = grid_.cells[face.neighbor_id_];
        const auto ass_face = face.GetNeighborAssociatedFace(grid_);
        const auto& adj_face = adj_cell.faces_[ass_face];

        auto& cur_face_ori = cell_face_orientations_[cell.local_id_][f];

        const double adj_mu = omega.Dot(adj_face.normal_);
        if (adj_mu > tolerance) orientation = FOOUTGOING;
        else if (adj_mu < tolerance)
          orientation = FOINCOMING;

        switch (orientation)
        {
          case FOPARALLEL:
            cur_face_ori = FOPARALLEL;
            break;
          case FOINCOMING:
            cur_face_ori = FOOUTGOING;
            break;
          case FOOUTGOING:
            cur_face_ori = FOINCOMING;
            break;
        }
      } // if not face owned locally at all

      ++f;
    } // for face
  }

  //============================================= Make directed connections
  for (auto& cell : grid_.local_cells)
  {
    const uint64_t c = cell.local_id_;
    size_t f = 0;
    for (auto& face : cell.faces_)
    {
      const double mu = omega.Dot(face.normal_);
      //======================================= If outgoing determine if
      //                                        it is to a local cell
      if (cell_face_orientations_[cell.local_id_][f] == FOOUTGOING)
      {
        //================================ If it is a cell and not bndry
        if (face.has_neighbor_)
        {
          //========================= If it is in the current location
          if (face.IsNeighborLocal(grid_))
          {
            double weight = mu * face.ComputeFaceArea(grid_);
            cell_successors[c].insert(std::make_pair(face.GetNeighborLocalID(grid_), weight));
          }
          else
            location_successors.insert(face.GetNeighborPartitionID(grid_));
        }
      }
      //======================================= If not outgoing determine
      //                                        what it is dependent on
      else
      {
        //================================if it is a cell and not bndry
        if (face.has_neighbor_ and not face.IsNeighborLocal(grid_))
          location_dependencies.insert(face.GetNeighborPartitionID(grid_));
      }
      ++f;
    } // for face
  }   // for cell
}

void
chi_mesh::sweep_management::SPDS::PrintedGhostedGraph() const
{
  constexpr double tolerance = 1.0e-16;

  for (int p = 0; p < Chi::mpi.process_count; ++p)
  {
    if (p == Chi::mpi.location_id)
    {
      std::cout << "// Printing directed graph for location " << p << ":\n";
      std::cout << "digraph DG {\n";
      std::cout << "  splines=\"FALSE\";\n";
      std::cout << "  rankdir=\"LR\";\n\n";
      std::cout << "  /* Vertices */\n";

      for (const auto& cell : grid_.local_cells)
      {
        std::cout << "  " << cell.global_id_ << " [shape=\"circle\"]\n";

        for (const auto& face : cell.faces_)
        {
          if (face.has_neighbor_ and (not grid_.IsCellLocal(face.neighbor_id_)))
            std::cout << "  " << face.neighbor_id_
                      << " [shape=\"circle\", style=filled fillcolor=red "
                         "fontcolor=white] //proc="
                      << grid_.cells[face.neighbor_id_].partition_id_ << "\n";
        }
      }

      std::cout << "\n"
                << "  /* Edges */\n";
      for (const auto& cell : grid_.local_cells)
      {
        for (const auto& face : cell.faces_)
        {
          if (face.has_neighbor_ and (cell.global_id_ > face.neighbor_id_))
          {
            if (omega_.Dot(face.normal_) > tolerance)
              std::cout << "  " << cell.global_id_ << " -> " << face.neighbor_id_ << "\n";
            else if (omega_.Dot(face.normal_) < tolerance)
              std::cout << "  " << face.neighbor_id_ << " -> " << cell.global_id_ << "\n";
          } // if outgoing
        }
      }
      std::cout << "\n";
      const auto ghost_ids = grid_.cells.GetGhostGlobalIDs();
      for (const uint64_t global_id : ghost_ids)
      {
        const auto& cell = grid_.cells[global_id];
        for (const auto& face : cell.faces_)
        {
          if (face.has_neighbor_ and (cell.global_id_ > face.neighbor_id_) and
              grid_.IsCellLocal(face.neighbor_id_))
          {
            if (omega_.Dot(face.normal_) > tolerance)
              std::cout << "  " << cell.global_id_ << " -> " << face.neighbor_id_ << "\n";
            else if (omega_.Dot(face.normal_) < tolerance)
              std::cout << "  " << face.neighbor_id_ << " -> " << cell.global_id_ << "\n";
          } // if outgoing
        }
      }
      std::cout << "}\n";

    } // if current location
    Chi::mpi.Barrier();
  } // for p
}
