#include "diffusion_solver.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Initialization of common to all solver types.*/
void chi_diffusion::Solver::InitializeCommonItems()
{
  const std::string fname = "chi_diffusion::Solver::InitializeCommonItems";
  grid_ptr_ = chi_mesh::GetCurrentHandler().GetGrid();

  if (grid_ptr_ == nullptr)
    throw std::logic_error(fname + " No grid defined.");

  auto globl_unique_bndry_ids = grid_ptr_->GetDomainUniqueBoundaryIDs();

  const auto& grid_boundary_id_map = grid_ptr_->GetBoundaryIDMap();
  for (uint64_t bndry_id : globl_unique_bndry_ids)
  {
    if (grid_boundary_id_map.count(bndry_id) == 0)
      throw std::logic_error(fname + ": Boundary id " +
      std::to_string(bndry_id) + " does not have a name-assignment.");

    const auto& bndry_name = grid_boundary_id_map.at(bndry_id);
    if (boundary_preferences_.find(bndry_name) != boundary_preferences_.end())
    {
      BoundaryInfo bndry_info = boundary_preferences_.at(bndry_name);
      auto& bndry_vals = bndry_info.second;
      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting:
        {
          boundaries_.insert(
            std::make_pair(bndry_id, new
                           chi_diffusion::BoundaryReflecting()));
          Chi::log.Log() << "Boundary \"" << bndry_name << "\" set to reflecting.";
          break;
        }
        case BoundaryType::Dirichlet:
        {
          if (bndry_vals.empty()) bndry_vals.resize(1,0.0);
          boundaries_.insert(
            std::make_pair(bndry_id, new
            chi_diffusion::BoundaryDirichlet(bndry_vals[0])));
          Chi::log.Log() << "Boundary \"" << bndry_name << "\" set to dirichlet.";
          break;
        }
        case BoundaryType::Robin:
        {
          if (bndry_vals.size()<3) bndry_vals.resize(3,0.0);
          boundaries_.insert(
            std::make_pair(bndry_id, new
                           chi_diffusion::BoundaryRobin(bndry_vals[0],
                                                        bndry_vals[1],
                                                        bndry_vals[2])));
          Chi::log.Log() << "Boundary \"" << bndry_name << "\" set to robin.";
          break;
        }
        case BoundaryType::Vacuum:
        {
          boundaries_.insert(
            std::make_pair(bndry_id, new
                           chi_diffusion::BoundaryRobin(0.25, 0.5, 0.0)));
          Chi::log.Log() << "Boundary \"" << bndry_name << "\" set to vacuum.";
          break;
        }
        case BoundaryType::Neumann:
        {
          if (bndry_vals.size()<3) bndry_vals.resize(3,0.0);
          boundaries_.insert(
            std::make_pair(bndry_id, new
            chi_diffusion::BoundaryRobin(bndry_vals[0],
                                         bndry_vals[1],
                                         bndry_vals[2])));
          Chi::log.Log() << "Boundary \"" << bndry_name << "\" set to neumann.";
          break;
        }
      }//switch boundary type
    }
    else
    {
      boundaries_.insert(
        std::make_pair(bndry_id, new
                       chi_diffusion::BoundaryDirichlet()));
      Chi::log.Log0Verbose1()
        << "No boundary preference found for boundary \"" << bndry_name
        << "\" Dirichlet boundary added with zero boundary value.";
    }
  }//for neighbor_id_

  common_items_initialized_ = true;
}


