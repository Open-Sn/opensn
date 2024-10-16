// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

int
SPDS::MapLocJToPrelocI(int locJ) const
{
  CALI_CXX_MARK_SCOPE("SPDS::MapLocJToPrelocI");

  for (int i = 0; i < location_dependencies_.size(); i++)
  {
    if (location_dependencies_[i] == locJ)
    {
      return i;
    }
  }

  for (int i = 0; i < delayed_location_dependencies_.size(); i++)
  {
    if (delayed_location_dependencies_[i] == locJ)
    {
      return -(i + 1);
    }
  }

  log.LogAllError() << "SPDS Invalid mapping encountered in MapLocJToPrelocI.";
  Exit(EXIT_FAILURE);
  return 0;
}

int
SPDS::MapLocJToDeplocI(int locJ) const
{
  CALI_CXX_MARK_SCOPE("SPDS::MapLocJToDeploc");

  for (int i = 0; i < location_successors_.size(); i++)
  {
    if (location_successors_[i] == locJ)
    {
      return i;
    }
  }

  log.LogAllError() << "SPDS Invalid mapping encountered in MapLocJToDeplocI.";
  Exit(EXIT_FAILURE);
  return 0;
}

void
SPDS::PopulateCellRelationships(const Vector3& omega,
                                std::set<int>& location_dependencies,
                                std::set<int>& location_successors,
                                std::vector<std::set<std::pair<int, double>>>& cell_successors)
{
  CALI_CXX_MARK_SCOPE("SPDS::PopulateCellRelationships");

  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations_.assign(grid_.local_cells.size(), {});
  for (auto& cell : grid_.local_cells)
    cell_face_orientations_[cell.local_id].assign(cell.faces.size(), FOPARALLEL);

  for (auto& cell : grid_.local_cells)
  {
    size_t f = 0;
    for (auto& face : cell.faces)
    {
      // Determine if the face is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal);

      bool owns_face = true;
      if (face.has_neighbor and cell.global_id > face.neighbor_id)
        owns_face = false;

      if (owns_face)
      {
        if (mu > tolerance)
          orientation = FOOUTGOING;
        else if (mu < tolerance)
          orientation = FOINCOMING;

        cell_face_orientations_[cell.local_id][f] = orientation;

        if (face.has_neighbor and grid_.IsCellLocal(face.neighbor_id))
        {
          const auto& adj_cell = grid_.cells[face.neighbor_id];
          const auto ass_face = face.GetNeighborAssociatedFace(grid_);
          auto& adj_face_ori = cell_face_orientations_[adj_cell.local_id][ass_face];

          switch (orientation)
          {
            case FOPARALLEL:
              adj_face_ori = FOPARALLEL;
              break;
            case FOINCOMING:
              adj_face_ori = FOOUTGOING;
              break;
            case FOOUTGOING:
              adj_face_ori = FOINCOMING;
              break;
          }
        }
      } // if face owned
      else if (face.has_neighbor and not grid_.IsCellLocal(face.neighbor_id))
      {
        const auto& adj_cell = grid_.cells[face.neighbor_id];
        const auto ass_face = face.GetNeighborAssociatedFace(grid_);
        const auto& adj_face = adj_cell.faces[ass_face];

        auto& cur_face_ori = cell_face_orientations_[cell.local_id][f];

        const double adj_mu = omega.Dot(adj_face.normal);
        if (adj_mu > tolerance)
          orientation = FOOUTGOING;
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

  // Make directed connections
  for (auto& cell : grid_.local_cells)
  {
    const uint64_t c = cell.local_id;
    size_t f = 0;
    for (auto& face : cell.faces)
    {
      const double mu = omega.Dot(face.normal);
      // If outgoing determine if it is to a local cell
      if (cell_face_orientations_[cell.local_id][f] == FOOUTGOING)
      {
        // If it is a cell and not bndry
        if (face.has_neighbor)
        {
          // If it is in the current location
          if (face.IsNeighborLocal(grid_))
          {
            double weight = mu * face.ComputeFaceArea(grid_);
            cell_successors[c].insert(std::make_pair(face.GetNeighborLocalID(grid_), weight));
          }
          else
            location_successors.insert(face.GetNeighborPartitionID(grid_));
        }
      }
      // If not outgoing determine what it is dependent on
      else
      {
        // if it is a cell and not bndry
        if (face.has_neighbor and not face.IsNeighborLocal(grid_))
          location_dependencies.insert(face.GetNeighborPartitionID(grid_));
      }
      ++f;
    } // for face
  }   // for cell
}

void
SPDS::PrintedGhostedGraph() const
{
  constexpr double tolerance = 1.0e-16;

  for (int p = 0; p < opensn::mpi_comm.size(); ++p)
  {
    if (p == opensn::mpi_comm.rank())
    {
      std::cout << "// Printing directed graph for location " << p << ":\n";
      std::cout << "digraph DG {\n";
      std::cout << "  splines=\"FALSE\";\n";
      std::cout << "  rankdir=\"LR\";\n\n";
      std::cout << "  /* Vertices */\n";

      for (const auto& cell : grid_.local_cells)
      {
        std::cout << "  " << cell.global_id << " [shape=\"circle\"]\n";

        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (not grid_.IsCellLocal(face.neighbor_id)))
            std::cout << "  " << face.neighbor_id
                      << " [shape=\"circle\", style=filled fillcolor=red "
                         "fontcolor=white] //proc="
                      << grid_.cells[face.neighbor_id].partition_id << "\n";
        }
      }

      std::cout << "\n"
                << "  /* Edges */\n";
      for (const auto& cell : grid_.local_cells)
      {
        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (cell.global_id > face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell.global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < tolerance)
              std::cout << "  " << face.neighbor_id << " -> " << cell.global_id << "\n";
          } // if outgoing
        }
      }
      std::cout << "\n";
      const auto ghost_ids = grid_.cells.GetGhostGlobalIDs();
      for (const uint64_t global_id : ghost_ids)
      {
        const auto& cell = grid_.cells[global_id];
        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (cell.global_id > face.neighbor_id) and
              grid_.IsCellLocal(face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell.global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < tolerance)
              std::cout << "  " << face.neighbor_id << " -> " << cell.global_id << "\n";
          } // if outgoing
        }
      }
      std::cout << "}\n";

    } // if current location
    opensn::mpi_comm.barrier();
  } // for p
}

} // namespace opensn
