#include "lbs_solver.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Sets a value to the zeroth (scalar) moment of the vector.*/
void
lbs::LBSSolver::SetPhiVectorScalarValues(std::vector<double>& phi_vector, double value)
{
  const size_t first_grp = groups_.front().id_;
  const size_t final_grp = groups_.back().id_;

  const auto& sdm = *discretization_;

  typedef const int64_t cint64_t;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_map = sdm.MapDOFLocal(cell,
                                         i,
                                         flux_moments_uk_man_,
                                         /*m*/ 0,
                                         /*g*/ 0);

      double* phi = &phi_vector[dof_map];

      for (size_t g = first_grp; g <= final_grp; ++g)
        phi[g] = value;
    } // for node i
  }   // for cell
}

//###################################################################
/**Scales a flux moment vector. For sweep methods the delayed angular
 * fluxes will also be scaled.*/
void
lbs::LBSSolver::ScalePhiVector(PhiSTLOption which_phi, double value)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  chi_math::Scale(*y_ptr, value);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void
lbs::LBSSolver::SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x, &x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          x_ref[index] = (*y_ptr)[mapping + g]; // Offset on purpose
        }                                       // for g
      }                                         // for moment
    }                                           // for dof
  }                                             // for cell

  VecRestoreArray(x, &x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void
lbs::LBSSolver::SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset,
                                                  Vec x_src,
                                                  PhiSTLOption which_phi)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromGSPETScVec");
  }

  const double* x_ref;
  VecGetArrayRead(x_src, &x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          (*y_ptr)[mapping + g] = x_ref[index];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell

  VecRestoreArrayRead(x_src, &x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void
lbs::LBSSolver::GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                              const std::vector<double>& x_src,
                                              std::vector<double>& y)
{
  int gsi = groupset.groups_.front().id_;
  size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          y[mapping + g] = x_src[mapping + g];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void
lbs::LBSSolver::GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                              PhiSTLOption from_which_phi,
                                              PhiSTLOption to_which_phi)
{
  std::vector<double>* y_ptr;
  switch (to_which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  std::vector<double>* x_src_ptr;
  switch (from_which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      x_src_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      x_src_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  int gsi = groupset.groups_.front().id_;
  size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          (*y_ptr)[mapping + g] = (*x_src_ptr)[mapping + g];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell
}

//###################################################################
/**Assembles a vector for a given group span from a source vector.*/
void
lbs::LBSSolver::SetGroupScopedPETScVecFromPrimarySTLvector(int first_group_id,
                                                           int last_group_id,
                                                           Vec x,
                                                           const std::vector<double>& y)
{
  double* x_ref;
  VecGetArray(x, &x_ref);

  int gsi = first_group_id;
  int gsf = last_group_id;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          x_ref[index] = y[mapping + g]; // Offset on purpose
        }                                // for g
      }                                  // for moment
    }                                    // for dof
  }                                      // for cell

  VecRestoreArray(x, &x_ref);
}

//###################################################################
/**Assembles a vector for a given groupset from a source vector.*/
void
lbs::LBSSolver::SetPrimarySTLvectorFromGroupScopedPETScVec(int first_group_id,
                                                           int last_group_id,
                                                           Vec x_src,
                                                           std::vector<double>& y)
{
  const double* x_ref;
  VecGetArrayRead(x_src, &x_ref);

  int gsi = first_group_id;
  int gsf = last_group_id;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          y[mapping + g] = x_ref[index];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell

  VecRestoreArrayRead(x_src, &x_ref);
}

//###################################################################
/**Assembles a PETSc vector from multiple groupsets.*/
void
lbs::LBSSolver::SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& gs_ids,
                                                       Vec x,
                                                       PhiSTLOption which_phi)
{
  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetMultiGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x, &x_ref);

  int64_t index = -1;
  for (int gs_id : gs_ids)
  {
    const auto& groupset = groupsets_.at(gs_id);

    int gsi = groupset.groups_.front().id_;
    int gsf = groupset.groups_.back().id_;
    int gss = gsf - gsi + 1;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      auto& transport_view = cell_transport_views_[cell.local_id_];

      for (int i = 0; i < cell.vertex_ids_.size(); i++)
      {
        for (int m = 0; m < num_moments_; m++)
        {
          size_t mapping = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; g++)
          {
            index++;
            x_ref[index] = (*y_ptr)[mapping + g]; // Offset on purpose
          }                                       // for g
        }                                         // for moment
      }                                           // for dof
    }                                             // for cell
  }                                               // for groupset id

  VecRestoreArray(x, &x_ref);
}

//###################################################################
/**Disassembles a multiple Groupset PETSc vector STL vectors.*/
void
lbs::LBSSolver::SetPrimarySTLvectorFromMultiGSPETScVecFrom(const std::vector<int>& gs_ids,
                                                           Vec x_src,
                                                           PhiSTLOption which_phi)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromMultiGSPETScVecFrom");
  }

  const double* x_ref;
  VecGetArrayRead(x_src, &x_ref);

  int64_t index = -1;
  for (int gs_id : gs_ids)
  {
    const auto& groupset = groupsets_.at(gs_id);

    int gsi = groupset.groups_.front().id_;
    int gsf = groupset.groups_.back().id_;
    int gss = gsf - gsi + 1;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      auto& transport_view = cell_transport_views_[cell.local_id_];

      for (int i = 0; i < cell.vertex_ids_.size(); i++)
      {
        for (int m = 0; m < num_moments_; m++)
        {
          size_t mapping = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; g++)
          {
            index++;
            (*y_ptr)[mapping + g] = x_ref[index];
          } // for g
        }   // for moment
      }     // for dof
    }       // for cell
  }         // for groupset id

  VecRestoreArrayRead(x_src, &x_ref);
}