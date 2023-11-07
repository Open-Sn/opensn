#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/sweep_utilities/sweep_namespace.h"
#include <fstream>

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/cell/cell.h"

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"

#include "framework/mesh/sweep_utilities/spds/spds.h"

void
chi_mesh::sweep_management::PrintSweepOrdering(chi_mesh::sweep_management::SPDS* sweep_order,
                                               MeshContinuumPtr vol_continuum)
{
  //  double polar = sweep_order->polar;
  //  double azimuthal = sweep_order->azimuthal;
  //
  //  chi_mesh::MeshHandler* cur_handler = GetCurrentHandler();
  //  chi_mesh::VolumeMesher* mesher     = cur_handler->volume_mesher;
  //
  //  // Making
  //  containers std::vector<int> rank_of_cell; for (int c=0;
  //  c<vol_continuum->cells.size(); c++)
  //  {
  //    rank_of_cell.push_back(0);
  //  }
  //
  //  chi_mesh::Vector3 omega;
  //  omega.x = sin(polar)*cos(azimuthal);
  //  omega.y = sin(polar)*sin(azimuthal);
  //  omega.z = cos(polar);
  //
  //
  //  // Traverse the
  //  graph int max_rank=-1;
  //
  //  for (int ci=0; ci<sweep_order->spls_->item_id.size(); ci++)
  //  {
  //    int cell_index = sweep_order->spls_->item_id[ci];
  //    chi_mesh::Cell* cell =  vol_continuum->cells[cell_index];
  //
  //    for (int e=0; e < cell->faces.size(); e++)
  //    {
  //      // Determine if the face
  //      //                                        is incident
  //      bool is_incoming = false;
  //      double dot_normal = omega.Dot(cell->faces[e].normal);
  //      if (dot_normal<0.0) {is_incoming = true;}
  //
  //      // If incoming determine if
  //      //                                        it is locally dependent
  //      if (is_incoming)
  //      {
  //        int adj_index = cell->faces[e].neighbor;
  //
  //        for (int lc=0;
  //             lc<vol_continuum->local_cell_glob_indices.size(); lc++)
  //        {
  //          if (adj_index == vol_continuum->local_cell_glob_indices[lc])
  //          {
  //            if (rank_of_cell[cell_index]==0)
  //            {
  //              rank_of_cell[cell_index] = rank_of_cell[adj_index]+1;
  //              printf("Sweep cell %d, rank
  //              %d\n",ci,rank_of_cell[cell_index]); if
  //              (max_rank<rank_of_cell[cell_index])
  //              {
  //                max_rank = rank_of_cell[cell_index];
  //              }
  //            }
  //            break;
  //          }
  //        }
  //      }
  //
  //    }//for face
  //  }
  //
  //  printf("Max rank=%d\n",max_rank);
  //
  //  // Sort ranks into
  //  //                                                        flags
  //  std::vector<int>* ranked_cells;
  //  std::vector<std::vector<int>*> ranks;
  //
  //  for (int i=0; i<(max_rank+1); i++)
  //  {
  //    ranked_cells = new std::vector<int>;
  //    ranks.push_back(ranked_cells);
  //  }
  //
  //  for (int c=0; c<vol_continuum->cells.size(); c++)
  //  {
  //    int rank = rank_of_cell[c];
  //    ranks[rank]->push_back(c);
  //  }
  //
  //  //for (int sp=0; sp<sweep_order->spls_.size(); sp++)
  //  for (int sp=0; sp<(max_rank+1); sp++)
  //  //for (int sp=0; sp<1; sp++)
  //  {
  //    std::string file_name("SweepMesh");
  //    file_name = file_name + std::to_string(sp)+std::string(".py");
  //    vol_continuum->ExportCellsToPython(
  //      file_name.c_str(),true,
  //      ranks[sp],1);
  //  }
}
