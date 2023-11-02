#include "FiniteVolume.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/SpatialDiscretization/CellMappings/FiniteVolumeMapping.h"

namespace chi_math::spatial_discretization
{

void
FiniteVolume::CreateCellMappings()
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization_FV::"
                                     "CreateCellMappings";

  auto MakeCellMapping = [this, fname](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
      case chi_mesh::CellType::POLYGON:
      case chi_mesh::CellType::POLYHEDRON:
      {
        typedef std::vector<std::vector<int>> FaceDofMapping;
        mapping = make_unique<cell_mapping::FiniteVolumeMapping>(
          ref_grid_, cell, cell.centroid_, FaceDofMapping(cell.faces_.size(), {-1}));
        break;
      }
      default:
        throw std::logic_error(std::string(fname) +
                               std::string(": Invalid cell type encountered."));
    }
    return mapping;
  };

  for (const auto& cell : ref_grid_.local_cells)
    cell_mappings_.push_back(MakeCellMapping(cell));

  const auto ghost_ids = ref_grid_.cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(ref_grid_.cells[ghost_id]);
    nb_cell_mappings_.insert(std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}

} // namespace chi_math::spatial_discretization
