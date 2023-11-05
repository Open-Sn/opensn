#include "opensn/framework/math/SpatialDiscretization/SpatialDiscretization.h"

#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "opensn/framework/math/PETScUtils/petsc_utils.h"

#include "opensn/framework/logging/chi_log.h"

namespace chi_math
{

SpatialDiscretization::SpatialDiscretization(const chi_mesh::MeshContinuum& grid,
                                             CoordinateSystemType cs_type,
                                             SDMType sdm_type)
  : UNITARY_UNKNOWN_MANAGER({std::make_pair(chi_math::UnknownType::SCALAR, 0)}),
    ref_grid_(grid),
    coord_sys_type_(cs_type),
    type_(sdm_type)
{
}

const CellMapping&
SpatialDiscretization::GetCellMapping(const chi_mesh::Cell& cell) const
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization::"
                                     "GetCellMapping";
  try
  {
    if (Grid().IsCellLocal(cell.global_id_)) return *cell_mappings_.at(cell.local_id_);
    else
      return *nb_cell_mappings_.at(cell.global_id_);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(std::string(fname) + ": Failed to obtain cell mapping.");
  }
}

SpatialDiscretizationType
SpatialDiscretization::Type() const
{
  return type_;
}

const chi_mesh::MeshContinuum&
SpatialDiscretization::Grid() const
{
  return ref_grid_;
}

CoordinateSystemType
SpatialDiscretization::GetCoordinateSystemType() const
{
  return coord_sys_type_;
}

size_t
SpatialDiscretization::GetNumLocalDOFs(const UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return local_base_block_size_ * N;
}

size_t
SpatialDiscretization::GetNumGlobalDOFs(const UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return globl_base_block_size_ * N;
}

size_t
SpatialDiscretization::GetNumLocalAndGhostDOFs(const UnknownManager& unknown_manager) const
{
  return GetNumLocalDOFs(unknown_manager) + GetNumGhostDOFs(unknown_manager);
}

size_t
SpatialDiscretization::GetCellNumNodes(const chi_mesh::Cell& cell) const
{
  return GetCellMapping(cell).NumNodes();
}

const std::vector<chi_mesh::Vector3>&
SpatialDiscretization::GetCellNodeLocations(const chi_mesh::Cell& cell) const
{
  return GetCellMapping(cell).GetNodeLocations();
}

std::pair<std::set<uint32_t>, std::set<uint32_t>>
SpatialDiscretization::MakeCellInternalAndBndryNodeIDs(const chi_mesh::Cell& cell) const
{
  const auto& cell_mapping = GetCellMapping(cell);
  const size_t num_faces = cell.faces_.size();
  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== Determine which nodes are on the
  //                                       boundary
  std::set<uint32_t> boundary_nodes;
  for (size_t f = 0; f < num_faces; ++f)
  {
    if (not cell.faces_[f].has_neighbor_)
    {
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
        boundary_nodes.insert(cell_mapping.MapFaceNode(f, fi));
    }
  } // for f

  //====================================== Determine non-boundary nodes
  std::set<uint32_t> internal_nodes;
  for (size_t i = 0; i < num_nodes; ++i)
    if (boundary_nodes.find(i) == boundary_nodes.end()) internal_nodes.insert(i);

  return {internal_nodes, boundary_nodes};
}

// ###################################################################
/**For each cell, for each face of that cell, for each node on that face,
 * maps to which local node on the adjacent cell that node position corresponds.
 *
\param tolerance double. Tolerance to use to determine if two node locations are
                         equal. [Default: 1.0e-12]
 *
 * For example consider two adjacent quadrilaterals.

\verbatim
o--------o--------o                       o--------o
|3      2|3      2|                       |    2   |
|  101   |   102  | , face ids for both:  |3      1|
|0      1|0      1|                       |    0   |
o--------o--------o                       o--------o
internal face for cell 101 is face-1, ccw orientated
--o
 1|
  |
 0|
--o
internal face for cell 102 is face-3, ccw orientated
o-
|0
|
|1
o-

mapping[101][1][0] = 0
mapping[101][1][1] = 3

mapping[102][3][0] = 2
mapping[102][3][1] = 1
\endverbatim

*/
std::vector<std::vector<std::vector<int>>>
chi_math::SpatialDiscretization::MakeInternalFaceNodeMappings(
  const double tolerance /*=1.0e-12*/) const
{
  typedef std::vector<int> FaceAdjMapping;
  typedef std::vector<FaceAdjMapping> PerFaceAdjMapping;
  typedef std::vector<PerFaceAdjMapping> CellAdjMapping;

  const auto& grid = this->ref_grid_;

  CellAdjMapping cell_adj_mapping;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = this->GetCellMapping(cell);
    const auto& node_locations = cell_mapping.GetNodeLocations();
    const size_t num_faces = cell.faces_.size();

    PerFaceAdjMapping per_face_adj_mapping;

    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      const auto num_face_nodes = cell_mapping.NumFaceNodes(f);
      FaceAdjMapping face_adj_mapping(num_face_nodes, -1);
      if (face.has_neighbor_)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id_];
        const auto& adj_cell_mapping = this->GetCellMapping(adj_cell);
        const auto& adj_node_locations = adj_cell_mapping.GetNodeLocations();
        const size_t adj_num_nodes = adj_cell_mapping.NumNodes();

        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const auto& ivec3 = node_locations[i];

          for (size_t ai = 0; ai < adj_num_nodes; ++ai)
          {
            const auto& aivec3 = adj_node_locations[ai];
            if ((ivec3 - aivec3).NormSquare() < tolerance)
            {
              face_adj_mapping[fi] = static_cast<int>(ai);
              break;
            }
          } // for ai
          if (face_adj_mapping[fi] < 0) throw std::logic_error("Face node mapping failed");
        } // for fi
      }   // if internal face

      per_face_adj_mapping.push_back(std::move(face_adj_mapping));
    } // for face

    cell_adj_mapping.push_back(std::move(per_face_adj_mapping));
  } // for cell

  return cell_adj_mapping;
}

// ###################################################################
/**Copy part of vector A to vector B. Suppose vector A's entries are
 * managed `chi_math::UnknownManager` A (`uk_manA`) and that the
 * entries of the vector B are managed by `chi_math::UnknownManager` B
 * (`uk_manB`). This function copies the entries associated with an unknown with
 * id `uk_id_A` in `uk_manA` from vector A to vector B such that the entries
 * in vector B are aligned with the entries of an unknown with id `uk_id_B` in
 * `uk_manB`. All the components are copied.
 *
\param from_vector Vector to copy from.
\param to_vector Vector to copy to.
\param from_vec_uk_structure Unknown manager for vector A.
\param from_vec_uk_id Unknown-id in unknown manager A.
\param to_vec_uk_structure Unknown manager for vector B.
\param to_vec_uk_id Unknown-id in unknown manager B.
 */
void
SpatialDiscretization::CopyVectorWithUnknownScope(const std::vector<double>& from_vector,
                                                  std::vector<double>& to_vector,
                                                  const UnknownManager& from_vec_uk_structure,
                                                  const unsigned int from_vec_uk_id,
                                                  const UnknownManager& to_vec_uk_structure,
                                                  const unsigned int to_vec_uk_id) const
{
  const std::string fname = "chi_math::SpatialDiscretization::"
                            "CopyVectorWithUnknownScope";
  const auto& ukmanF = from_vec_uk_structure;
  const auto& ukmanT = to_vec_uk_structure;
  const auto& ukidF = from_vec_uk_id;
  const auto& ukidT = to_vec_uk_id;
  try
  {
    const auto& ukA = from_vec_uk_structure.unknowns_.at(from_vec_uk_id);
    const auto& ukB = to_vec_uk_structure.unknowns_.at(to_vec_uk_id);

    if (ukA.num_components_ != ukB.num_components_)
      throw std::logic_error(fname + " Unknowns do not have the "
                                     "same number of components");

    const size_t num_comps = ukA.num_components_;

    for (const auto& cell : ref_grid_.local_cells)
    {
      const auto& cell_mapping = this->GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        for (size_t c = 0; c < num_comps; ++c)
        {
          const int64_t fmap = MapDOFLocal(cell, i, ukmanF, ukidF, c);
          const int64_t imap = MapDOFLocal(cell, i, ukmanT, ukidT, c);

          to_vector[imap] = from_vector[fmap];
        } // for component c
      }   // for node i
    }     // for cell
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(fname + ": either from_vec_uk_id or to_vec_uk_id is "
                                    "out of range for its respective "
                                    "unknown manager.");
  }
}

// ###################################################################
/**Develops a localized view of a petsc vector.*/
void
SpatialDiscretization::LocalizePETScVector(Vec petsc_vector,
                                           std::vector<double>& local_vector,
                                           const chi_math::UnknownManager& unknown_manager) const
{
  size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);

  chi_math::PETScUtils::CopyVecToSTLvector(petsc_vector, local_vector, num_local_dofs);
}

// ###################################################################
/**Develops a localized view of a petsc vector.*/
void
SpatialDiscretization::LocalizePETScVectorWithGhosts(
  Vec petsc_vector,
  std::vector<double>& local_vector,
  const chi_math::UnknownManager& unknown_manager) const
{
  size_t num_local_dofs = GetNumLocalAndGhostDOFs(unknown_manager);

  chi_math::PETScUtils::CopyVecToSTLvectorWithGhosts(petsc_vector, local_vector, num_local_dofs);
}

double
SpatialDiscretization::CartesianSpatialWeightFunction(const chi_mesh::Vector3& point)
{
  return 1.0;
}

double
SpatialDiscretization::CylindricalRZSpatialWeightFunction(const chi_mesh::Vector3& point)
{
  return 2.0 * M_PI * point[0];
}

double
SpatialDiscretization::Spherical1DSpatialWeightFunction(const chi_mesh::Vector3& point)
{
  const double r = point[2];
  return 4.0 * M_PI * r * r;
}

SpatialDiscretization::SpatialWeightFunction
SpatialDiscretization::GetSpatialWeightingFunction() const
{
  switch (coord_sys_type_)
  {
    case CoordinateSystemType::CARTESIAN:
      return CartesianSpatialWeightFunction;
    case CoordinateSystemType::CYLINDRICAL:
      return CylindricalRZSpatialWeightFunction;
    case CoordinateSystemType::SPHERICAL:
      return Spherical1DSpatialWeightFunction;
    case CoordinateSystemType::UNDEFINED:
    default:
      ChiLogicalError("Coordinate system undefined.");
  }
}

} // namespace chi_math
