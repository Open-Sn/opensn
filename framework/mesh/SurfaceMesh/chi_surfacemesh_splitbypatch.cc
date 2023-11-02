#include "chi_surfacemesh.h"

/**Splits the surface by patch.*/
void
chi_mesh::SurfaceMesh::SplitByPatch(std::vector<chi_mesh::SurfaceMesh*>& patches)
{
  typedef std::vector<chi_mesh::Face> FaceList;
  typedef std::vector<FaceList*> FaceListCollection;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Copy all faces from surface
  FaceList unsorted_faces;
  FaceList::iterator cur_face;
  for (cur_face = this->faces_.begin(); cur_face != this->faces_.end(); cur_face++)
  {
    unsorted_faces.push_back((*cur_face));
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize co-planar
  // collection
  FaceListCollection co_planar_lists;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Seed the first collection
  FaceList* first_list = new FaceList;
  co_planar_lists.push_back(first_list);
  first_list->push_back(unsorted_faces.back());
  unsorted_faces.pop_back();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reverse iterate over unsorted
  FaceList::reverse_iterator us_face;
  for (us_face = unsorted_faces.rbegin(); us_face != unsorted_faces.rend(); us_face++)
  {
    bool matchFound = false;

    //====================================== Check if it can be matched
    FaceListCollection::iterator existing_list;
    for (existing_list = co_planar_lists.begin(); existing_list != co_planar_lists.end();
         existing_list++)
    {
      for (cur_face = (*existing_list)->begin(); cur_face != (*existing_list)->end(); cur_face++)
      {
        chi_mesh::Normal n1 = cur_face->geometric_normal;
        chi_mesh::Normal n2 = us_face->geometric_normal;

        if (fabs(n1.Dot(n2)) > (1.0 - 1.0e-4))
        {
          (*existing_list)->push_back(unsorted_faces.back());
          unsorted_faces.pop_back();
          matchFound = true;
          break;
        }
      }
      if (matchFound) { break; }
    }

    //====================================== If no match found, create new list
    if (!matchFound)
    {
      printf("New list created\n");
      FaceList* new_list = new FaceList;
      new_list->push_back(unsorted_faces.back());
      unsorted_faces.pop_back();
      co_planar_lists.push_back(new_list);
    }
  }

  printf("Number of co-planar collections=%lu\n", co_planar_lists.size());

  FaceListCollection patch_list_collection;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over co_planar lists
  FaceListCollection::iterator existing_list;
  for (existing_list = co_planar_lists.begin(); existing_list != co_planar_lists.end();
       existing_list++)
  {
    //================================= Inner patch collection
    FaceListCollection inner_patch_list_collection;

    //================================= Add all the faces of this collection
    //                                  to an unused list
    FaceList unused_faces;
    for (cur_face = (*existing_list)->begin(); cur_face != (*existing_list)->end(); cur_face++)
    {
      unused_faces.push_back((*cur_face));
    }

    //================================= Seed the first patch list
    FaceList* first_patch_list = new FaceList;
    inner_patch_list_collection.push_back(first_patch_list);
    first_patch_list->push_back(unused_faces.back());
    unused_faces.pop_back();

    //================================= Loop over unused faces
    while (unused_faces.size() > 0)
    {
      // printf("Unused faces=%d\n",unused_faces.size());
      bool updateMade = false;
      FaceList::iterator us_face_f; // Forward iterator
      for (us_face_f = unused_faces.begin(); us_face_f != unused_faces.end(); us_face_f++)
      {
        //============================= Try to to find a home
        FaceListCollection::iterator inner_list;
        for (inner_list = inner_patch_list_collection.begin();
             inner_list != inner_patch_list_collection.end();
             inner_list++)
        {
          for (cur_face = (*inner_list)->begin(); cur_face != (*inner_list)->end(); cur_face++)
          {
            //==================== Check if any vertices match
            for (int e = 0; e < 3; e++)
            {
              int vi = (*us_face_f).v_index[e];
              for (int e2 = 0; e2 < 3; e2++)
              {
                int vf = (*cur_face).v_index[e2];

                if (vf == vi)
                {
                  (*inner_list)->push_back(*us_face_f);
                  unused_faces.erase(us_face_f);
                  updateMade = true;
                  break;
                }
              }
              if (updateMade) { break; }
            }
            if (updateMade) { break; }
          }
          if (updateMade) { break; }
        }

        if (updateMade) { break; }
      }

      if (!updateMade)
      {
        FaceList* new_patch_list = new FaceList;
        inner_patch_list_collection.push_back(new_patch_list);
        new_patch_list->push_back(unused_faces.back());
        unused_faces.pop_back();
      }
    }

    //================================= Add inner patch lists to outer
    FaceListCollection::iterator inner_list;
    for (inner_list = inner_patch_list_collection.begin();
         inner_list != inner_patch_list_collection.end();
         inner_list++)
    {
      patch_list_collection.push_back(*inner_list);
    }
  }

  printf("Number of patches = %lu\n", patch_list_collection.size());

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create surfaces for each patch
  FaceListCollection::iterator outer_patch;
  for (outer_patch = patch_list_collection.begin(); outer_patch != patch_list_collection.end();
       outer_patch++)
  {

    chi_mesh::SurfaceMesh* new_surface = new chi_mesh::SurfaceMesh;

    std::vector<int*> vertex_mapping;

    for (cur_face = (*outer_patch)->begin(); cur_face != (*outer_patch)->end(); cur_face++)
    {
      //==================================== Copy the face
      chi_mesh::Face newFace = (*cur_face);

      //==================================== Copy and map vertices
      for (int e = 0; e < 3; e++)
      {
        int vi = newFace.v_index[e];

        //============================= Check if vertex already there
        bool already_there = false;
        int* already_there_mapping;
        std::vector<int*>::iterator cur_map;
        for (cur_map = vertex_mapping.begin(); cur_map != vertex_mapping.end(); cur_map++)
        {
          if ((*cur_map)[0] == vi)
          {
            already_there = true;
            already_there_mapping = (*cur_map);
            break;
          }
        }

        if (already_there)
        {
          //=========================== Update vertex index
          newFace.v_index[e] = already_there_mapping[1];
        }
        else
        {
          //=========================== Copy vertex
          chi_mesh::Vertex v = this->vertices_.at(vi);
          int* newMapping = new int[2];
          newMapping[0] = vi;
          newMapping[1] = new_surface->vertices_.size();

          new_surface->vertices_.push_back(v);
          vertex_mapping.push_back(newMapping);

          newFace.v_index[e] = newMapping[1];
        }

      } // for e
      newFace.e_index[0][0] = newFace.v_index[0];
      newFace.e_index[0][1] = newFace.v_index[1];

      newFace.e_index[1][0] = newFace.v_index[1];
      newFace.e_index[1][1] = newFace.v_index[2];

      newFace.e_index[2][0] = newFace.v_index[2];
      newFace.e_index[2][1] = newFace.v_index[0];

      for (int e = 0; e < 3; e++)
      {
        newFace.e_index[e][2] = -1;
        newFace.e_index[e][3] = -1;
      }

      new_surface->faces_.push_back(newFace);
    }
    new_surface->UpdateInternalConnectivity();
    // std::cout << (*new_surface);
    patches.push_back(new_surface);

    //    std::string fileName = "ZZMesh";
    //    fileName = fileName + std::to_string((int)patches.size());
    //    fileName = fileName + ".obj";
    //
    //    std::cout <<fileName <<"\n";
    //
    //    new_surface->ExportToOBJFile(fileName.c_str());
  }
}