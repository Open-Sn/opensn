#include"chi_surfacemesh.h"
#include<iostream>

bool chi_mesh::SurfaceMesh::CheckNegativeSense(double x, double y, double z)
{
  chi_mesh::Vector3 xyz = chi_mesh::Vector3(x, y, z);

  //======================================================= Loop through each face
  std::vector<chi_mesh::Face>::iterator cur_face;
  for (cur_face = this->faces_.begin();
          cur_face != this->faces_.end(); cur_face++)
  {
    //=========================================== Get a vertex (first one)
    chi_mesh::Vertex p;
    try{
      p = this->vertices_.at(cur_face->v_index[0]);
    }
    catch(const std::out_of_range& o){
      std::cout << "Invalid vertex handle" << std::endl;
    }

    //=========================================== Calculate dot product
    chi_mesh::Vector3 p_xyz = xyz - p;
    double dprod = cur_face->assigned_normal.Dot(p_xyz);

    if (dprod<0.0)
    {
      return true;
    }
  }

  return false;
}
