#include "diffusion_solver.h"

#include "physics/PhysicsMaterial/chi_physicsmaterial.h"
#include "physics/PhysicsMaterial/material_property_scalarvalue.h"
#include "physics/PhysicsMaterial/MultiGroupXS/multigroup_xs.h"
#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"

//###################################################################
/**Gets material properties various sources.*/
void
chi_diffusion::Solver::GetMaterialProperties(const chi_mesh::Cell& cell,
                                             int cell_dofs,
                                             std::vector<double>& diffCoeff,
                                             std::vector<double>& sourceQ,
                                             std::vector<double>& sigmaa,
                                             int group,
                                             int moment)
{
  uint64_t cell_glob_index = cell.global_id_;
  bool cell_is_local = (cell.partition_id_ == Chi::mpi.location_id);
  uint64_t cell_local_id = cell.local_id_;
  int mat_id = cell.material_id_;

  if (mat_id < 0)
  {
    Chi::log.Log0Error() << "Cell encountered with no material id. ";
    Chi::Exit(EXIT_FAILURE);
  }

  if (mat_id >= Chi::material_stack.size())
  {
    Chi::log.Log0Error() << "Cell encountered with material id pointing to "
                            "non-existing material.";
    Chi::Exit(EXIT_FAILURE);
  }

  auto property_map_D = basic_options_("property_map_D").IntegerValue();
  auto property_map_q = basic_options_("property_map_q").IntegerValue();
  auto property_map_sigma = basic_options_("property_map_sigma").IntegerValue();

  auto material = Chi::GetStackItemPtr(Chi::material_stack, mat_id, __FUNCTION__);

  //====================================== Process material properties
  diffCoeff.resize(cell_dofs, 1.0);
  sourceQ.resize(cell_dofs, 0.0);
  sigmaa.resize(cell_dofs, 0.0);

  //####################################################### REGULAR MATERIAL
  if (material_mode_ == DIFFUSION_MATERIALS_REGULAR)
  {
    // We absolutely need the diffusion coefficient so process error
    if ((property_map_D < 0) || (property_map_D >= material->properties_.size()))
    {
      Chi::log.Log0Error() << "Solver diffusion coefficient mapped to property index "
                           << property_map_D << " is not a valid index for material \""
                           << material->name_ << "\" id " << mat_id;
      Chi::Exit(EXIT_FAILURE);
    }

    // For now, we can only support scalar values so lets check that
    if (std::dynamic_pointer_cast<chi_physics::ScalarValue>(material->properties_[property_map_D]))
    {
      diffCoeff.assign(cell_dofs, material->properties_[property_map_D]->GetScalarValue());
    }
    else
    {
      Chi::log.Log0Error() << "Solver diffusion coefficient mapped to property index "
                           << property_map_D << " is not a valid property type"
                           << " for material \"" << material->name_ << "\" id " << mat_id
                           << ". Currently SCALAR_VALUE and THERMAL_CONDUCTIVITY are the "
                           << "only supported types.";
      Chi::Exit(EXIT_FAILURE);
    }

    if ((property_map_q < material->properties_.size()) && (property_map_q >= 0))
    {
      if (std::dynamic_pointer_cast<chi_physics::ScalarValue>(
            material->properties_[property_map_q]))
      {
        sourceQ.assign(cell_dofs, material->properties_[property_map_q]->GetScalarValue());
      }
      else
      {
        Chi::log.Log0Error() << "Source value mapped to property index " << property_map_q
                             << " is not a valid property type"
                             << " for material \"" << material->name_ << "\" id " << mat_id
                             << ". Currently SCALAR_VALUE is the "
                             << "only supported type.";
        Chi::Exit(EXIT_FAILURE);
      }
    }

    if (not((property_map_sigma < 0) || (property_map_sigma >= material->properties_.size())))
    {
      sigmaa.assign(cell_dofs, material->properties_[property_map_sigma]->GetScalarValue());
    }
  } // regular

  //####################################################### TRANSPORT XS D
  //                                                        TRANSPORT XS SIGA
  //                                                        SCALAR       Q
  else if (material_mode_ == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR)
  {
    //====================================== Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p = 0; p < material->properties_.size(); p++)
    {
      if (std::dynamic_pointer_cast<chi_physics::MultiGroupXS>(material->properties_[p]))
      {
        auto xs = std::static_pointer_cast<chi_physics::MultiGroupXS>(material->properties_[p]);

        diffCoeff.assign(cell_dofs, xs->DiffusionCoefficient()[group]);
        sigmaa.assign(cell_dofs, xs->SigmaRemoval()[group]);
        transportxs_found = true;
      }
    } // for properties

    if (!transportxs_found)
    {
      Chi::log.LogAllError() << "Diffusion Solver: Material encountered with no tranport xs"
                                " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      Chi::Exit(EXIT_FAILURE);
    }

    //====================================== Setting Q
    if ((property_map_q < material->properties_.size()) && (property_map_q >= 0))
    {
      if (std::dynamic_pointer_cast<chi_physics::ScalarValue>(
            material->properties_[property_map_q]))
      {
        sourceQ.assign(cell_dofs, material->properties_[property_map_q]->GetScalarValue());
      }
      else
      {
        Chi::log.Log0Error() << "Source value mapped to property index " << property_map_q
                             << " is not a valid property type"
                             << " for material \"" << material->name_ << "\" id " << mat_id
                             << ". Currently SCALAR_VALUE is the "
                             << "only supported type.";
        Chi::Exit(EXIT_FAILURE);
      }
    }
  } // transport xs TTR
  else
  {
    Chi::log.Log0Error() << "Diffusion Solver: Invalid material mode.";
    Chi::Exit(EXIT_FAILURE);
  }
}

//###################################################################
/**Update the field functions with the latest data.*/
void
chi_diffusion::Solver::UpdateFieldFunctions()
{
  Chi::log.LogAll() << "Updating field functions" << std::endl;
  auto& ff = *field_functions_.front();
  const auto& OneDofPerNode = discretization_->UNITARY_UNKNOWN_MANAGER;

  std::vector<double> data_vector;
  discretization_->LocalizePETScVector(x_, data_vector, OneDofPerNode);

  ff.UpdateFieldVector(data_vector);
}
