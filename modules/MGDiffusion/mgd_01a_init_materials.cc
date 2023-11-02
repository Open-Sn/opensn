#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "mesh/MeshHandler/chi_meshhandler.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearContinuous.h"
#include "physics/PhysicsMaterial/chi_physicsmaterial.h"
#include "physics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"
#include "physics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include <algorithm>
//============================================= assemble matrix A
void
mg_diffusion::Solver::Initialize_Materials(std::set<int>& material_ids)
{
  Chi::log.Log0Verbose1() << "Initializing Materials";

  std::stringstream materials_list;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process materials found
  const size_t num_physics_mats = Chi::material_stack.size();
  bool first_material_read = true;

  for (const int& mat_id : material_ids)
  {
    auto current_material = Chi::GetStackItemPtr(Chi::material_stack, mat_id, __FUNCTION__);
    materials_list << "Material id " << mat_id;

    //====================================== Check valid ids
    if (mat_id < 0)
    {
      throw std::logic_error("MG-diff-InitMaterials: Cells encountered with no assigned material.");
    }
    if (static_cast<size_t>(mat_id) >= num_physics_mats)
    {
      throw std::logic_error("MG-diff-InitMaterials: Cells encountered with "
                             "material id that matches no material in physics material library.");
    }

    //====================================== Extract properties
    using MatProperty = chi_physics::PropertyType;
    bool found_transport_xs = false;
    for (const auto& property : current_material->properties_)
    {
      if (property->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs = std::static_pointer_cast<chi_physics::MultiGroupXS>(property);
        matid_to_xs_map[mat_id] = transp_xs;
        found_transport_xs = true;
        if (first_material_read) num_groups_ = transp_xs->NumGroups();

      } // transport xs
      if (property->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source = std::static_pointer_cast<chi_physics::IsotropicMultiGrpSource>(property);

        if (mg_source->source_value_g_.size() < num_groups_)
        {
          Chi::log.LogAllWarning()
            << "MG-Diff-InitMaterials: Isotropic Multigroup source specified "
               "in "
            << "material \"" << current_material->name_ << "\" has fewer "
            << "energy groups than called for in the simulation. "
            << "Source will be ignored.";
        }
        else { matid_to_src_map[mat_id] = mg_source; }
      } // P0 source
    }   // for property

    //====================================== Check valid property
    if (!found_transport_xs)
    {
      Chi::log.LogAllError() << "MG-Diff-InitMaterials: Found no transport cross-section property "
                                "for "
                             << "material \"" << current_material->name_ << "\".";
      Chi::Exit(EXIT_FAILURE);
    }
    //====================================== Check number of groups legal
    if (matid_to_xs_map[mat_id]->NumGroups() != num_groups_)
    {
      Chi::log.LogAllError() << "MG-Diff-InitMaterials: Found material \""
                             << current_material->name_ << "\" has "
                             << matid_to_xs_map[mat_id]->NumGroups() << " groups and "
                             << "the simulation has " << num_groups_ << " groups. The material "
                             << "must have the same number of groups.";
      Chi::Exit(EXIT_FAILURE);
    }

    //====================================== Check number of moments
    if (matid_to_xs_map[mat_id]->ScatteringOrder() > 1)
    {
      Chi::log.Log0Warning() << "MG-Diff-InitMaterials: Found material \""
                             << current_material->name_ << "\" has a scattering order of "
                             << matid_to_xs_map[mat_id]->ScatteringOrder() << " and"
                             << " the simulation has a scattering order of One (MG-Diff)"
                             << " The higher moments will therefore not be used.";
    }

    materials_list << " number of moments " << matid_to_xs_map[mat_id]->ScatteringOrder() + 1
                   << "\n";

    first_material_read = false;
  } // for material id

  Chi::log.Log() << "Materials Initialized:\n" << materials_list.str() << "\n";

  Chi::mpi.Barrier();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute last fast group
  // initialize last fast group
  Chi::log.Log() << "Computing last fast group.";
  unsigned int lfg = num_groups_;

  if (num_groups_ > 1)
  {
    // loop over all materials
    for (const auto& mat_id_xs : matid_to_xs_map)
    {
      // get the P0 transfer matrix
      const auto& S = mat_id_xs.second->TransferMatrix(0);
      // loop over all row of the transfer matrix
      const int G = static_cast<int>(num_groups_);
      for (int g = G - 1; g >= 0; --g)
      {
        for (const auto& [row_g, gp, sigma_sm] : S.Row(g))
        {
          if ((std::fabs(sigma_sm) > 1e-10) && (gp > row_g))
            lfg = std::min(lfg, static_cast<unsigned int>(row_g));
        }
      }
    } // loop on materials
  }   // if num_groups>1

  last_fast_group_ = lfg;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute two-grid params
  do_two_grid_ = basic_options_("do_two_grid").BoolValue();
  if ((lfg == num_groups_) && do_two_grid_)
  {
    Chi::log.Log0Error() << "Two-grid is not possible with no upscattering.";
    do_two_grid_ = false;
    Chi::Exit(EXIT_FAILURE);
  }
  if (do_two_grid_)
  {
    Chi::log.Log() << "Compute_TwoGrid_Params";
    Compute_TwoGrid_Params();
  }
}
