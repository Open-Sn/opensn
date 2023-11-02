#include "chi_surfacemesh.h"

//#########################################################
/**Default constructor.*/
chi_mesh::SurfaceMesh::SurfaceMesh()
{
}

//#########################################################
/**Default destructor.*/
chi_mesh::SurfaceMesh::~SurfaceMesh()
{
  for (auto poly_face : poly_faces_)
  {
    delete poly_face;
  }

  poly_faces_.clear();
}

//#########################################################
std::ostream&
operator<<(std::ostream& os, chi_mesh::SurfaceMesh& that)
{
  std::vector<chi_mesh::Face>::const_iterator curface;
  for (curface = that.GetTriangles().begin(); curface != that.GetTriangles().end(); curface++)
  {
    long index = std::distance(that.GetTriangles().begin(), curface);
    os << "Face " << index << " v:";
    os << curface->v_index[0] << "->";
    os << curface->v_index[1] << "->";
    os << curface->v_index[2] << "  ";

    os << "e:";
    for (int e = 0; e < 3; e++)
    {
      os << "[" << curface->e_index[e][0] << ",";
      os << curface->e_index[e][1] << ",";
      os << curface->e_index[e][2] << ",";
      os << curface->e_index[e][3] << "]";
    }
    os << std::endl;
  }

  return os;
}