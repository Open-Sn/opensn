#include "AAH_FLUDSCommonData.h"

#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_mesh::sweep_management
{

void AAH_FLUDSCommonData::NonLocalIncidentMapping(const chi_mesh::Cell& cell,
                                                  const SPDS& spds)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();

  //=================================================== Loop over faces
  //           INCIDENT                                 but process
  //                                                    only incident faces
  for (short f=0; f < cell.faces_.size(); f++)
  {
    const CellFace&  face = cell.faces_[f];
    const auto& orientation = spds.CellFaceOrientations()[cell.local_id_][f];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident face
    if (orientation == FaceOrientation::INCOMING)
    {
      if ((face.has_neighbor_) and (!face.IsNeighborLocal(grid)) )
      {
        //============================== Find prelocI
        int locJ = face.GetNeighborPartitionID(grid);
        int prelocI = spds.MapLocJToPrelocI(locJ);

        //###########################################################
        if (prelocI >= 0)
        {
          //============================== Find the cell in prelocI cell views
          int ass_cell = -1;
          for (int c=0; c<prelocI_cell_views[prelocI].size(); c++)
          {
            if (prelocI_cell_views[prelocI][c].first == face.neighbor_id_)
            {
              ass_cell = c;
              break;
            }
          }
          if (ass_cell<0)
          {
            Chi::log.LogAll()
              << "Required predecessor cell not located in call to"
              << " InitializeBetaElements. locJ=" << locJ
              << " prelocI=" << prelocI
              << " cell=" << face.neighbor_id_;
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Find associated face
          std::set<int> cfvids(face.vertex_ids_.begin(),
                               face.vertex_ids_.end());
          CompactCellView* adj_cell_view =
            &prelocI_cell_views[prelocI][ass_cell];
          int ass_face = -1, af = -1;
          for (auto& adj_face : adj_cell_view->second)
          {
            ++af;
            bool face_matches = true;

            std::set<int> afvids(adj_face.second.begin(),
                                 adj_face.second.end());

            if (cfvids != afvids) face_matches = false;

            if (face_matches){ass_face = af; break;}
          }
          if (ass_face<0)
          {
            Chi::log.LogAll()
              << "Associated face not found in call to InitializeBetaElements";
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Map dofs
          std::pair<int,std::vector<int>> dof_mapping;
          dof_mapping.first = adj_cell_view->second[ass_face].first;
          std::vector<uint64_t>* ass_face_verts =
            &adj_cell_view->second[ass_face].second;
          for (int fv=0; fv < face.vertex_ids_.size(); fv++)
          {
            bool match_found = false;
            for (int afv=0; afv<ass_face_verts->size(); afv++)
            {
              if (face.vertex_ids_[fv] ==
                  ass_face_verts->operator[](afv))
              {
                match_found = true;
                dof_mapping.second.push_back(afv);
                break;
              }
            }

            if (!match_found)
            {
              Chi::log.LogAll()
                << "Associated vertex not found in call to InitializeBetaElements";
              Chi::Exit(EXIT_FAILURE);
            }
          }

          //============================== Push back final face info
          std::pair<int,std::pair<int,std::vector<int>>> inc_face_prelocI_info;
          inc_face_prelocI_info.first = prelocI;
          inc_face_prelocI_info.second =
            std::pair<int,std::vector<int>>(dof_mapping);

          std::pair<int,std::pair<int,std::vector<int>>> empty_delayed_info;
          empty_delayed_info.first = prelocI;

          nonlocal_inc_face_prelocI_slot_dof.push_back(inc_face_prelocI_info);
          delayed_nonlocal_inc_face_prelocI_slot_dof.push_back(empty_delayed_info);
        }//If not delayed predecessor
        //###########################################################
        else
        {
          int delayed_preLocI = abs(prelocI)-1;
          //============================== Find the cell in prelocI cell views
          int ass_cell = -1;
          for (int c=0; c<delayed_prelocI_cell_views[delayed_preLocI].size(); c++)
          {
            if (delayed_prelocI_cell_views[delayed_preLocI][c].first == face.neighbor_id_)
            {
              ass_cell = c;
              break;
            }
          }
          if (ass_cell<0)
          {
            Chi::log.LogAll()
              << "Required predecessor cell not located in call to"
              << " InitializeBetaElements. locJ=" << locJ
              << " delayed prelocI=" << delayed_preLocI
              << " cell=" << face.neighbor_id_;
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Find associated face
          CompactCellView* adj_cell_view =
            &delayed_prelocI_cell_views[delayed_preLocI][ass_cell];
          int ass_face = -1;
          for (int af=0; af<adj_cell_view->second.size(); af++)
          {
            bool face_matches = true;
            for (int afv=0; afv<adj_cell_view->second[af].second.size(); afv++)
            {
              bool match_found = false;
              for (int fv=0; fv < face.vertex_ids_.size(); fv++)
              {
                if (adj_cell_view->second[af].second[afv] ==
                    face.vertex_ids_[fv])
                {
                  match_found = true;
                  break;
                }
              }

              if (!match_found){face_matches = false; break;}
            }

            if (face_matches){ass_face = af; break;}
          }
          if (ass_face<0)
          {
            Chi::log.LogAll()
              << "Associated face not found in call to InitializeBetaElements";
            Chi::Exit(EXIT_FAILURE);
          }

          //============================== Map dofs
          std::pair<int,std::vector<int>> dof_mapping;
          dof_mapping.first = adj_cell_view->second[ass_face].first;
          std::vector<uint64_t>* ass_face_verts =
            &adj_cell_view->second[ass_face].second;
          for (int fv=0; fv < face.vertex_ids_.size(); fv++)
          {
            bool match_found = false;
            for (int afv=0; afv<ass_face_verts->size(); afv++)
            {
              if (face.vertex_ids_[fv] ==
                  ass_face_verts->operator[](afv))
              {
                match_found = true;
                dof_mapping.second.push_back(afv);
                break;
              }
            }

            if (!match_found)
            {
              Chi::log.LogAll()
                << "Associated vertex not found in call to InitializeBetaElements";
              Chi::Exit(EXIT_FAILURE);
            }
          }

          //============================== Push back final face info
          std::pair<int,std::pair<int,std::vector<int>>> inc_face_prelocI_info;
          inc_face_prelocI_info.first = delayed_preLocI;
          inc_face_prelocI_info.second =
            std::pair<int,std::vector<int>>(dof_mapping);

          std::pair<int,std::pair<int,std::vector<int>>> delayed_info;
          delayed_info.first = -(delayed_preLocI + 1);

          delayed_nonlocal_inc_face_prelocI_slot_dof.push_back(inc_face_prelocI_info);
          nonlocal_inc_face_prelocI_slot_dof.push_back(delayed_info);
        }//If delayed predecessor

      }//if not local and not boundary
    }//if incident
  }//for incindent f
}

} // namespace chi_mesh::sweep_management