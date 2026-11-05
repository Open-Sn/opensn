// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <algorithm>
#include <vector>
#include <queue>
#include <stdexcept>

namespace opensn
{

SPDSFaceNeighborInfoVec
BuildSPDSFaceNeighborInfo(const Mesh& grid)
{
  SPDSFaceNeighborInfoVec info;
  info.resize(grid.GetLocalCellCount());
  for (const auto& cell : grid.GetLocalCells())
  {
    auto& cell_info = info[cell->local_id];
    cell_info.resize(cell->faces.size());
    for (size_t f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      auto& fi = cell_info[f];
      fi.has_neighbor = face.has_neighbor;
      fi.owns_face = true;
      if (face.has_neighbor and cell->global_id > face.neighbor_id)
        fi.owns_face = false;

      if (face.has_neighbor)
      {
        fi.neighbor_adj_face = face.GetNeighborAdjacentFaceIndex(&grid);
        if (face.IsNeighborLocal(&grid))
        {
          fi.neighbor_local = true;
          fi.neighbor_local_id = face.GetNeighborLocalID(&grid);
        }
        else
        {
          fi.neighbor_partition_id = face.GetNeighborPartitionID(&grid);
        }
      }
    }
  }

  return info;
}

int
SPDS::MapLocJToPrelocI(int locJ) const
{

  for (int i = 0; i < location_dependencies_.size(); ++i)
  {
    if (location_dependencies_[i] == locJ)
    {
      return i;
    }
  }

  for (int i = 0; i < delayed_location_dependencies_.size(); ++i)
  {
    if (delayed_location_dependencies_[i] == locJ)
    {
      return -(i + 1);
    }
  }

  throw std::runtime_error("SPDS: Invalid mapping encountered in MapLocJToPrelocI");
  return 0;
}

int
SPDS::MapLocJToDeplocI(int locJ) const
{

  for (int i = 0; i < location_successors_.size(); ++i)
  {
    if (location_successors_[i] == locJ)
    {
      return i;
    }
  }

  throw std::runtime_error("SPDS: Invalid mapping encountered in MapLocJToDeplocI");
  return 0;
}

void
SPDS::PopulateCellRelationships(
  const Vector3& omega,
  const SPDSFaceNeighborInfoVec& face_info,
  std::vector<int>& location_dependencies,
  std::vector<int>& location_successors,
  std::vector<std::vector<std::pair<std::uint32_t, double>>>& cell_successors)
{

  constexpr double tolerance = FACE_ORIENTATION_TOLERANCE;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations_.assign(grid_->GetLocalCellCount(), {});
  for (const auto& cell : grid_->GetLocalCells())
    cell_face_orientations_[cell->local_id].assign(cell->faces.size(), FOPARALLEL);

  for (const auto& cell : grid_->GetLocalCells())
  {
    size_t f = 0;
    for (auto& face : cell->faces)
    {
      // Determine if the face is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal);

      const auto& finfo = face_info[cell->local_id][f];
      const bool owns_face = finfo.owns_face;

      if (owns_face)
      {
        if (mu > tolerance)
          orientation = FOOUTGOING;
        else if (mu < -tolerance)
          orientation = FOINCOMING;

        cell_face_orientations_[cell->local_id][f] = orientation;

        if (finfo.has_neighbor and finfo.neighbor_local)
        {
          const auto& adj_cell = grid_->GetLocalCell(finfo.neighbor_local_id);
          const auto adj_face_idx = finfo.neighbor_adj_face;
          auto& adj_face_ori = cell_face_orientations_[adj_cell.local_id][adj_face_idx];

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
      else if (finfo.has_neighbor and not finfo.neighbor_local)
      {
        const auto& adj_cell = grid_->GetGlobalCell(face.neighbor_id);
        const auto adj_face_idx = finfo.neighbor_adj_face;
        const auto& adj_face = adj_cell.faces[adj_face_idx];

        auto& cur_face_ori = cell_face_orientations_[cell->local_id][f];

        const double adj_mu = omega.Dot(adj_face.normal);
        if (adj_mu > tolerance)
          orientation = FOOUTGOING;
        else if (adj_mu < -tolerance)
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
  for (const auto& cell : grid_->GetLocalCells())
  {
    size_t f = 0;
    for (auto& face : cell->faces)
    {
      const double mu = omega.Dot(face.normal);
      // If outgoing determine if it is to a local cell
      if (cell_face_orientations_[cell->local_id][f] == FOOUTGOING)
      {
        // If it is a cell and not bndry
        if (face_info[cell->local_id][f].has_neighbor)
        {
          // If it is in the current location
          if (face_info[cell->local_id][f].neighbor_local)
          {
            const auto weight = mu * face.area;
            cell_successors[cell->local_id].emplace_back(
              face_info[cell->local_id][f].neighbor_local_id, weight);
          }
          else
            location_successors.push_back(face_info[cell->local_id][f].neighbor_partition_id);
        }
      }
      // If not outgoing determine what it is dependent on
      else if (cell_face_orientations_[cell->local_id][f] == FOINCOMING)
      {
        // if it is a cell and not bndry
        if (face_info[cell->local_id][f].has_neighbor and
            not face_info[cell->local_id][f].neighbor_local)
          location_dependencies.push_back(face_info[cell->local_id][f].neighbor_partition_id);
      }
      ++f;
    } // for face
  } // for cell

  auto SortUnique = [](auto& values)
  {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
  };

  SortUnique(location_dependencies);
  SortUnique(location_successors);
  for (auto& successors : cell_successors)
    SortUnique(successors);
}

void
SPDS::PrintGhostedGraph() const
{
  constexpr double tolerance = FACE_ORIENTATION_TOLERANCE;

  for (int p = 0; p < opensn::mpi_comm.size(); ++p)
  {
    if (p == opensn::mpi_comm.rank())
    {
      std::cout << "// Printing directed graph for location " << p << ":\n";
      std::cout << "digraph DG {\n";
      std::cout << "  splines=\"FALSE\";\n";
      std::cout << "  rankdir=\"LR\";\n\n";
      std::cout << "  /* Vertices */\n";

      for (const auto& cell : grid_->GetLocalCells())
      {
        std::cout << "  " << cell->global_id << " [shape=\"circle\"]\n";

        for (const auto& face : cell->faces)
        {
          if (face.has_neighbor and (not grid_->IsCellLocal(face.neighbor_id)))
            std::cout << "  " << face.neighbor_id
                      << " [shape=\"circle\", style=filled fillcolor=red "
                         "fontcolor=white] //proc="
                      << grid_->GetGlobalCell(face.neighbor_id).partition_id << "\n";
        }
      }

      std::cout << "\n"
                << "  /* Edges */\n";
      for (const auto& cell : grid_->GetLocalCells())
      {
        for (const auto& face : cell->faces)
        {
          if (face.has_neighbor and (cell->global_id > face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell->global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < -tolerance)
              std::cout << "  " << face.neighbor_id << " -> " << cell->global_id << "\n";
          } // if outgoing
        }
      }
      std::cout << "\n";
      const auto ghost_ids = grid_->GetGhostGlobalIDs();
      for (const uint64_t global_id : ghost_ids)
      {
        const auto& cell = grid_->GetGlobalCell(global_id);
        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (cell.global_id > face.neighbor_id) and
              grid_->IsCellLocal(face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell.global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < -tolerance)
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
