#include "chi_volumemesher.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
chi_mesh::VolumeMesher::VolumeMesher(VolumeMesherType type) : type_(type)
{
}

//###################################################################
/** Sets the grid member of the volume mesher.*/
void
chi_mesh::VolumeMesher::SetContinuum(MeshContinuumPtr& grid)
{
  grid_ptr_ = grid;
}

//###################################################################
/** Gets the smart-pointer for the grid.*/
chi_mesh::MeshContinuumPtr&
chi_mesh::VolumeMesher::GetContinuum()
{
  return grid_ptr_;
}

//###################################################################
/**Sets grid attributes. This is normally a private member of the grid
 * but this class is a friend.*/
void
chi_mesh::VolumeMesher::SetGridAttributes(MeshAttributes new_attribs,
                                          std::array<size_t, 3> ortho_Nis /*={0,0,0}*/)
{
  grid_ptr_->SetAttributes(new_attribs, ortho_Nis);
}

//###################################################################
/**Gets the volume mesher's type.*/
chi_mesh::VolumeMesherType
chi_mesh::VolumeMesher::Type() const
{
  return type_;
}
