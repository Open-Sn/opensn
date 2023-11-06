#include "framework/mesh/SweepUtilities/AngleAggregation/angleaggregation.h"
#include "framework/mesh/SweepUtilities/SweepBoundary/boundary_reflecting.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/mpi/chi_mpi.h"

#define ExceptionReflectedAngleError                                                               \
  std::logic_error(fname + "Reflected angle not found for angle " + std::to_string(n) +            \
                   " with direction " + quadrature->omegas_[n].PrintStr() +                        \
                   ". This can happen for two reasons: i) A quadrature is used"                    \
                   " that is not symmetric about the axis associated with the "                    \
                   "reflected boundary, or ii) the reflecting boundary is not "                    \
                   "aligned with any reflecting axis of the quadrature.")

// ###################################################################
/** Sets up the angle-aggregation object. */
chi_mesh::sweep_management::AngleAggregation::AngleAggregation(
  const std::map<uint64_t, SweepBndryPtr>& in_sim_boundaries,
  size_t in_number_of_groups,
  size_t in_number_of_group_subsets,
  std::shared_ptr<chi_math::AngularQuadrature>& in_quadrature,
  chi_mesh::MeshContinuumPtr& in_grid)
{
  sim_boundaries = in_sim_boundaries;
  number_of_groups = in_number_of_groups;
  number_of_group_subsets = in_number_of_group_subsets;
  quadrature = in_quadrature;
  grid = in_grid;

  for (auto& bndry_id_cond : sim_boundaries)
    bndry_id_cond.second->Setup(*grid, *quadrature);

  is_setup = true;
}

// ###################################################################
/** Resets all the outgoing intra-location and inter-location
 * cyclic interfaces.*/
void
chi_mesh::sweep_management::AngleAggregation::ZeroOutgoingDelayedPsi()
{
  for (auto& angsetgrp : angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      for (auto& delayed_data : angset->GetFLUDS().DelayedPrelocIOutgoingPsi())
        chi_math::Set(delayed_data, 0.0);

  for (auto& angsetgrp : angle_set_groups)
    for (auto& angset : angsetgrp.AngleSets())
      chi_math::Set(angset->GetFLUDS().DelayedLocalPsi(), 0.0);
}

// ###################################################################
/** Resets all the incoming intra-location and inter-location
 * cyclic interfaces.*/
void
chi_mesh::sweep_management::AngleAggregation::ZeroIncomingDelayedPsi()
{
  //======================================== Opposing reflecting bndries
  for (const auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = 0.0;

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      chi_math::Set(angle_set->GetFLUDS().DelayedLocalPsiOld(), 0.0);

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        chi_math::Set(loc_vector, 0.0);
}

// ###################################################################
/** Initializes reflecting boundary conditions. */
void
chi_mesh::sweep_management::AngleAggregation::InitializeReflectingBCs()
{
  const std::string fname = "chi_mesh::sweep_management::AngleAggregation";
  const double epsilon = 1.0e-8;

  bool reflecting_bcs_initialized = false;

  const chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  const chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);

  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      size_t tot_num_angles = quadrature->abscissae_.size();
      size_t num_local_cells = grid->local_cells.size();
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      const auto& normal = rbndry.Normal();

      rbndry.GetReflectedAngleIndexMap().resize(tot_num_angles, -1);
      rbndry.GetAngleReadyFlags().resize(tot_num_angles,
                                         std::vector<bool>(number_of_group_subsets, false));

      //========================================= Determine reflected angle
      //                                          and check that it is within
      //                                          the quadrature
      typedef chi_mesh::Vector3 Vec3;
      for (int n = 0; n < tot_num_angles; ++n)
      {
        const Vec3& omega_n = quadrature->omegas_[n];
        Vec3 omega_reflected;

        switch (rbndry.CoordType())
        {
          case chi_math::CoordinateSystemType::SPHERICAL:
            omega_reflected = -1.0 * omega_n;
            break;
          case chi_math::CoordinateSystemType::CYLINDRICAL:
          {
            // left, top and bottom is regular reflecting
            if (std::fabs(normal.Dot(jhat)) > 0.999999 or normal.Dot(ihat) < -0.999999)
              omega_reflected = omega_n - 2.0 * normal * omega_n.Dot(normal);
            // right derive their normal from omega_n
            else if (normal.Dot(ihat) > 0.999999)
            {
              Vec3 normal_star;
              if (omega_n.Dot(normal) > 0.0)
                normal_star = Vec3(omega_n.x, 0.0, omega_n.z).Normalized();
              else
                normal_star = Vec3(-omega_n.x, 0.0, -omega_n.y).Normalized();

              omega_reflected = omega_n - 2.0 * normal_star * omega_n.Dot(normal_star);
            }
          }
          break;
          case chi_math::CoordinateSystemType::CARTESIAN:
          default:
            omega_reflected = omega_n - 2.0 * normal * omega_n.Dot(normal);
            break;
        }

        auto& index_map = rbndry.GetReflectedAngleIndexMap();
        for (int nstar = 0; nstar < tot_num_angles; ++nstar)
          if (omega_reflected.Dot(quadrature->omegas_[nstar]) > (1.0 - epsilon))
          {
            index_map[n] = nstar;
            break;
          }

        if (index_map[n] < 0) throw ExceptionReflectedAngleError;
      }

      //========================================= Initialize storage for all
      //                                          outbound directions
      auto& heteroflux_new = rbndry.GetHeteroBoundaryFluxNew();
      auto& heteroflux_old = rbndry.GetHeteroBoundaryFluxOld();
      heteroflux_new.clear();
      heteroflux_old.clear();
      heteroflux_new.resize(tot_num_angles);
      for (int n = 0; n < tot_num_angles; ++n)
      {
        // Only continue if omega is outgoing
        if (quadrature->omegas_[n].Dot(rbndry.Normal()) < 0.0) continue;

        //================================== For cells
        auto& cell_vec = heteroflux_new[n];
        cell_vec.resize(num_local_cells);
        for (const auto& cell : grid->local_cells)
        {
          const uint64_t c = cell.local_id_;

          //=========================== Check cell on ref bndry
          bool on_ref_bndry = false;
          for (const auto& face : cell.faces_)
          {
            if ((not face.has_neighbor_) and (face.normal_.Dot(rbndry.Normal()) > 0.999999))
            {
              on_ref_bndry = true;
              break;
            }
          }
          if (not on_ref_bndry) continue;

          //=========================== If cell on ref bndry
          cell_vec[c].resize(cell.faces_.size());
          int f = 0;
          for (const auto& face : cell.faces_)
          {
            if ((not face.has_neighbor_) and (face.normal_.Dot(rbndry.Normal()) > 0.999999))
            {
              cell_vec[c][f].clear();
              cell_vec[c][f].resize(face.vertex_ids_.size(),
                                    std::vector<double>(number_of_groups, 0.0));
            }
            ++f;
          }
        } // for cells
      }   // for angles

      //========================================= Determine if boundary is
      //                                          opposing reflecting
      // The boundary with the smallest bid will
      // be marked as "opposing-reflecting" while
      // the other one will be just a regular
      // reflecting boundary
      for (const auto& [otherbid, otherbndry] : sim_boundaries)
      {
        if (bid == otherbid) continue;
        if (not otherbndry->IsReflecting()) continue;

        const auto& otherRbndry = dynamic_cast<const BoundaryReflecting&>(*otherbndry);

        if (rbndry.Normal().Dot(otherRbndry.Normal()) < (0.0 - epsilon))
          if (bid < otherbid) rbndry.SetOpposingReflected(true);
      }

      if (rbndry.IsOpposingReflected())
        rbndry.GetHeteroBoundaryFluxOld() = rbndry.GetHeteroBoundaryFluxNew();

      reflecting_bcs_initialized = true;
    } // if reflecting
  }   // for bndry

  if (reflecting_bcs_initialized)
    Chi::log.Log0Verbose1() << "Reflecting boundary conditions initialized.";
}

// ###################################################################
/** Returns a pair of numbers containing the number of
 * delayed angular unknowns both locally and globally, respectively. */
std::pair<size_t, size_t>
chi_mesh::sweep_management::AngleAggregation::GetNumDelayedAngularDOFs()
{
  //======================================== Check if this is already developed
  if (num_ang_unknowns_avail) return number_angular_unknowns;

  //======================================== If not developed
  size_t local_ang_unknowns = 0;

  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                local_ang_unknowns += dofvec.size();

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      local_ang_unknowns += angle_set->GetFLUDS().DelayedLocalPsi().size();

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
        local_ang_unknowns += loc_vector.size();

  size_t global_ang_unknowns = 0;
  MPI_Allreduce(
    &local_ang_unknowns, &global_ang_unknowns, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chi::mpi.comm);

  number_angular_unknowns = {local_ang_unknowns, global_ang_unknowns};

  num_ang_unknowns_avail = true;
  return number_angular_unknowns;
}

// ###################################################################
/** Assembles angular unknowns into the reference vector. */
void
chi_mesh::sweep_management::AngleAggregation::AppendNewDelayedAngularDOFsToArray(int64_t& index,
                                                                                 double* x_ref)
{
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                {
                  index++;
                  x_ref[index] = val;
                }

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsi())
      {
        index++;
        x_ref[index] = val;
      }

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
        for (auto val : loc_vector)
        {
          index++;
          x_ref[index] = val;
        }
}

// ###################################################################
/** Assembles angular unknowns into the reference vector. */
void
chi_mesh::sweep_management::AngleAggregation::AppendOldDelayedAngularDOFsToArray(int64_t& index,
                                                                                 double* x_ref)
{
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                {
                  index++;
                  x_ref[index] = val;
                }

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsiOld())
      {
        index++;
        x_ref[index] = val;
      }

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        for (auto val : loc_vector)
        {
          index++;
          x_ref[index] = val;
        }
}

// ###################################################################
/** Assembles angular unknowns into the reference vector. */
void
chi_mesh::sweep_management::AngleAggregation::SetOldDelayedAngularDOFsFromArray(int64_t& index,
                                                                                const double* x_ref)
{
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                {
                  index++;
                  val = x_ref[index];
                }

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsiOld())
      {
        index++;
        val = x_ref[index];
      }

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        for (auto& val : loc_vector)
        {
          index++;
          val = x_ref[index];
        }
}

// ###################################################################
/** Assembles angular unknowns into the reference vector. */
void
chi_mesh::sweep_management::AngleAggregation::SetNewDelayedAngularDOFsFromArray(int64_t& index,
                                                                                const double* x_ref)
{
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                {
                  index++;
                  val = x_ref[index];
                }

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsi())
      {
        index++;
        val = x_ref[index];
      }

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
        for (auto& val : loc_vector)
        {
          index++;
          val = x_ref[index];
        }
}

// ###################################################################
/**Gets the current values of the angular unknowns as an STL vector.*/
std::vector<double>
chi_mesh::sweep_management::AngleAggregation::GetNewDelayedAngularDOFsAsSTLVector()
{
  std::vector<double> psi_vector;

  auto psi_size = GetNumDelayedAngularDOFs();
  psi_vector.reserve(psi_size.first);

  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                  psi_vector.push_back(val);

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsi())
        psi_vector.push_back(val);

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
        for (auto val : loc_vector)
          psi_vector.push_back(val);

  return psi_vector;
}

// ###################################################################
/**Gets the current values of the angular unknowns as an STL vector.*/
void
chi_mesh::sweep_management::AngleAggregation::SetNewDelayedAngularDOFsFromSTLVector(
  const std::vector<double>& stl_vector)
{
  auto psi_size = GetNumDelayedAngularDOFs();
  size_t stl_size = stl_vector.size();
  if (stl_size != psi_size.first)
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": STL-vector size "
                           "is incompatible with number angular unknowns stored "
                           "in the angle-aggregation object.");

  size_t index = 0;
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxNew())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = stl_vector[index++];

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsi())
        val = stl_vector[index++];

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi())
        for (auto& val : loc_vector)
          val = stl_vector[index++];
}

// ###################################################################
/**Gets the current values of the angular unknowns as an STL vector.*/
std::vector<double>
chi_mesh::sweep_management::AngleAggregation::GetOldDelayedAngularDOFsAsSTLVector()
{
  std::vector<double> psi_vector;

  auto psi_size = GetNumDelayedAngularDOFs();
  psi_vector.reserve(psi_size.first);

  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto val : dofvec)
                  psi_vector.push_back(val);

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto val : angle_set->GetFLUDS().DelayedLocalPsiOld())
        psi_vector.push_back(val);

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        for (auto val : loc_vector)
          psi_vector.push_back(val);

  return psi_vector;
}

// ###################################################################
/**Gets the current values of the angular unknowns as an STL vector.*/
void
chi_mesh::sweep_management::AngleAggregation::SetOldDelayedAngularDOFsFromSTLVector(
  const std::vector<double>& stl_vector)
{
  auto psi_size = GetNumDelayedAngularDOFs();
  size_t stl_size = stl_vector.size();
  if (stl_size != psi_size.first)
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": STL-vector size "
                           "is incompatible with number angular unknowns stored "
                           "in the angle-aggregation object.");

  size_t index = 0;
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        for (auto& angle : rbndry.GetHeteroBoundaryFluxOld())
          for (auto& cellvec : angle)
            for (auto& facevec : cellvec)
              for (auto& dofvec : facevec)
                for (auto& val : dofvec)
                  val = stl_vector[index++];

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& val : angle_set->GetFLUDS().DelayedLocalPsiOld())
        val = stl_vector[index++];

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      for (auto& loc_vector : angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld())
        for (auto& val : loc_vector)
          val = stl_vector[index++];
}

// ###################################################################
/**Copies the old delayed angular fluxes to the new.*/
void
chi_mesh::sweep_management::AngleAggregation::SetDelayedPsiOld2New()
{
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        rbndry.GetHeteroBoundaryFluxNew() = rbndry.GetHeteroBoundaryFluxOld();

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedLocalPsi() = angle_set->GetFLUDS().DelayedLocalPsiOld();

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi() =
        angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld();
}

// ###################################################################
/**Copies the new delayed angular fluxes to the old.*/
void
chi_mesh::sweep_management::AngleAggregation::SetDelayedPsiNew2Old()
{
  //======================================== Opposing reflecting bndries
  for (auto& [bid, bndry] : sim_boundaries)
  {
    if (bndry->IsReflecting())
    {
      auto& rbndry = (BoundaryReflecting&)(*bndry);

      if (rbndry.IsOpposingReflected())
        rbndry.GetHeteroBoundaryFluxOld() = rbndry.GetHeteroBoundaryFluxNew();

    } // if reflecting
  }   // for bndry

  //======================================== Intra-cell cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedLocalPsiOld() = angle_set->GetFLUDS().DelayedLocalPsi();

  //======================================== Inter location cycles
  for (auto& as_group : angle_set_groups)
    for (auto& angle_set : as_group.AngleSets())
      angle_set->GetFLUDS().DelayedPrelocIOutgoingPsiOld() =
        angle_set->GetFLUDS().DelayedPrelocIOutgoingPsi();
}
