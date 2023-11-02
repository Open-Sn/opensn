#ifndef _chi_meshface_h
#define _chi_meshface_h

//######################################################### Struct
/**Data structure for a triangular face.*/
struct chi_mesh::Face
{
  int v_index[3];
  int n_index[3];
  int vt_index[3];
  int e_index[3][4];

  chi_mesh::Normal geometric_normal;
  chi_mesh::Normal assigned_normal;
  chi_mesh::Vertex face_centroid;

  bool invalidated;

  Face()
  {
    for (int k=0;k<3;k++)
    {
      v_index[k]=-1;
      n_index[k]=-1;
      vt_index[k]=-1;
      e_index[k][0]=-1;
      e_index[k][1]=-1;
      e_index[k][2]=-1;
      e_index[k][3]=-1;
      invalidated = false;
    }
  }

  void SetIndices(int a, int b, int c)
  {
    v_index[0] = a;
    v_index[1] = b;
    v_index[2] = c;

    e_index[0][0] = a; e_index[0][1] = b;
    e_index[1][0] = b; e_index[1][1] = c;
    e_index[2][0] = c; e_index[2][1] = a;
  }

  Face& operator=(const Face& that)
  {
    for (int k=0;k<3;k++)
    {
      v_index[k] =that.v_index[k] ;
      n_index[k] =that.n_index[k] ;
      vt_index[k]=that.vt_index[k];
      e_index[k][0]=that.e_index[k][0];
      e_index[k][1]=that.e_index[k][1];
      e_index[k][2]=that.e_index[k][2];
      e_index[k][3]=that.e_index[k][3];
    }
    geometric_normal = that.geometric_normal;
    assigned_normal  = that.assigned_normal;
    invalidated = that.invalidated;
    return *this;
  }

};

#define NEIGHBOR 0

//######################################################### Struct
/**Data structure for a polygon face.

 edges\n
 An array of 4 integers.\n
 [0] = Vertex index of edge start.\n
 [1] = Vertex index of edge end.\n
 [2] = Index of the face adjoining this edge (not the current face).
       -1 if not connected to anything,-1*boundary_index if connected
       to a boundary.\n
 [3] = Edge number of adjoining face. -1 if not connected
       to anything. 0 if a boundary.\n
   \n
   \n
  face_indices\n
  [0] = Index of the adjoining cell. -1 if not connected to anything.
        -1*boundary_index if connected to a boundary.\n
  [1] = Face number of adjoining cell. -1 if not connected
       to anything. 0 if a boundary.\n
  [2] = Partition ID of adjecent cell.\n

 */
struct chi_mesh::PolyFace
{
  std::vector<int> v_indices;
  std::vector<int*> edges;
  int face_indices[3];

  chi_mesh::Normal geometric_normal;
  chi_mesh::Vertex face_centroid;

  bool invalidated;

  PolyFace()
  {
    invalidated = false;
  }

  ~PolyFace()
  {
    for (auto edge : edges) delete [] edge;
  }

};


#endif
