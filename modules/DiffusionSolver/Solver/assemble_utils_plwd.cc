#include "diffusion_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Still searching for a reference for this.
 *
 * For Polygons:
 * Defined from paper  \n
 * Turcksin B, Ragusa J, "Discontinuous diffusion synthetic acceleration
 * for S_n transport on 2D arbitrary polygonal meshes", Journal of
 * Computational Physics 274, pg 356-369, 2014.\n
 * \n
 * Nv = Number of vertices. If Nv <= 4 then the perimeter parameter
 * should be replaced by edge length.*/
double
chi_diffusion::Solver::HPerpendicular(const chi_mesh::Cell& cell,
                                      const UnitIntegralContainer& fe_intgrl_values,
                                      unsigned int f)
{
  double hp = 1.0;

  size_t Nf = cell.faces_.size();
  size_t Nv = cell.vertex_ids_.size();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    const auto& v0 = grid_ptr_->vertices[cell.vertex_ids_[0]];
    const auto& v1 = grid_ptr_->vertices[cell.vertex_ids_[1]];

    hp = (v1 - v0).Norm() / 2.0;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    //    Nv = 4;
    const chi_mesh::CellFace& face = cell.faces_[f];

    uint64_t v0i = face.vertex_ids_[0];
    uint64_t v1i = face.vertex_ids_[1];

    const auto& v0 = grid_ptr_->vertices[v0i];
    const auto& v1 = grid_ptr_->vertices[v1i];

    double perimeter = (v1 - v0).Norm();

    double area = 0.0;
    for (int i = 0; i < fe_intgrl_values.NumNodes(); i++)
      area += fe_intgrl_values.IntV_shapeI(i);

    hp = area / perimeter;

    //    if (Nv == 3)
    //      hp = 2*area/perimeter;
    //    else if (Nv == 4)
    //      hp = area/perimeter;
    //    else //Nv > 4
    //    {
    //      if (Nv%2 == 0)
    //        hp = 4*area/perimeter;
    //      else
    //      {
    //        hp = 2*area/perimeter;
    //        hp += sqrt(2*area / Nv*sin(2*M_PI/Nv));
    //      }
    //    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    double volume = 0.0;
    for (int i = 0; i < fe_intgrl_values.NumNodes(); i++)
      volume += fe_intgrl_values.IntV_shapeI(i);

    double area = 0.0;
    for (int fr = 0; fr < Nf; fr++)
      for (int i = 0; i < Nv; i++)
        area += fe_intgrl_values.IntS_shapeI(fr, i);

    if (Nf == 4) // Tet
      hp = 3 * volume / area;
    else if (Nf == 6 && Nv == 8) // Hex
      hp = volume / area;
    else // Polyhedron
      hp = 6 * volume / area;
  } // Polyhedron
  else
  {
    Chi::log.LogAllError() << "Unsupported cell type in call to HPerpendicular";
    Chi::Exit(EXIT_FAILURE);
  }

  return hp;
}

/**Given a global node index, returns the local cell-node it's associated on the
 * referenced cell.*/
uint64_t
chi_diffusion::Solver::MapCellLocalNodeIDFromGlobalID(const chi_mesh::Cell& cell,
                                                      uint64_t node_global_id)
{
  size_t imap = 0;
  bool map_found = false;
  for (size_t ai = 0; ai < cell.vertex_ids_.size(); ai++)
  {
    if (node_global_id == cell.vertex_ids_[ai])
    {
      imap = ai;
      map_found = true;
      break;
    }
  }

  if (not map_found) throw std::logic_error(std::string(__FUNCTION__) + ": Mapping failure.");

  return imap;
}

/**Given the face index on the current cell, finds the
 * corresponding face index on the adjacent cell.*/
unsigned int
chi_diffusion::Solver::MapCellFace(const chi_mesh::Cell& cur_cell,
                                   const chi_mesh::Cell& adj_cell,
                                   unsigned int f)
{
  const auto& ccface = cur_cell.faces_[f]; // current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids_)
    ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af = 0; af < adj_cell.faces_.size(); af++)
  {
    const auto& acface = adj_cell.faces_[af]; // adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids_)
      acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  } // for adj faces

  if (not map_found) throw std::logic_error(std::string(__FUNCTION__) + ": Mapping failure.");

  return (unsigned int)fmap;
}
