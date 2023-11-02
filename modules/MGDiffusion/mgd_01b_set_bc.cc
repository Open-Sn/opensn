#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "mg_diffusion_bndry.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"

#include "physics/PhysicsMaterial/chi_physicsmaterial.h"

//============================================= assemble matrix A
void
mg_diffusion::Solver::Set_BCs(const std::vector<uint64_t>& globl_unique_bndry_ids)
{
  Chi::log.Log0Verbose1() << "Setting Boundary Conditions";

  uint64_t max_boundary_id = 0;
  for (const auto& id : globl_unique_bndry_ids)
    max_boundary_id = std::max(id, max_boundary_id);

  Chi::log.Log() << "Max boundary id identified: " << max_boundary_id;

  for (int bndry = 0; bndry < (max_boundary_id + 1); ++bndry)
  {
    if (boundary_preferences_.find(bndry) != boundary_preferences_.end())
    {
      BoundaryInfo bndry_info = boundary_preferences_.at(bndry);
      auto& bndry_vals = bndry_info.second;

      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting: // ------------- REFLECTING
        {
          boundaries_.push_back({BoundaryType::Reflecting});
          Chi::log.Log() << "Boundary " << bndry << " set to reflecting.";
          break;
        }
        case BoundaryType::Robin: // ------------- ROBIN
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Robin needs 3 values in bndry vals.");
          boundaries_.push_back(Boundary{BoundaryType::Robin, bndry_vals});
          Chi::log.Log() << "Boundary " << bndry << " set to robin.";
          break;
        }
        case BoundaryType::Vacuum: // ------------- VACUUM
        {
          auto ng = mg_diffusion::Solver::num_groups_;
          std::vector<double> a_values(ng, 0.25);
          std::vector<double> b_values(ng, 0.5);
          std::vector<double> f_values(ng, 0.0);
          boundaries_.push_back(Boundary{BoundaryType::Robin, {a_values, b_values, f_values}});
          Chi::log.Log() << "Boundary " << bndry << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann: // ------------- NEUMANN
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Neumann needs 3 values in bndry vals.");
          boundaries_.push_back(Boundary{BoundaryType::Robin, bndry_vals});
          Chi::log.Log() << "Boundary " << bndry << " set to neumann.";
          break;
        }
      } // switch boundary type
    }
    else
    {
      auto ng = mg_diffusion::Solver::num_groups_;
      std::vector<double> a_values(ng, 0.25);
      std::vector<double> b_values(ng, 0.5);
      std::vector<double> f_values(ng, 0.0);
      boundaries_.push_back({BoundaryType::Robin, {a_values, b_values, f_values}});
      Chi::log.Log0Verbose1() << "No boundary preference found for boundary index " << bndry
                              << "Vacuum boundary added as default.";
    }
  } // for bndry
}
