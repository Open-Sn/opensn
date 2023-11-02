#include "lbs_solver.h"

#include "Tools/lbs_bndry_func_lua.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define mk_shrd(x) std::make_shared<x>
#define SweepVaccuumBndry chi_mesh::sweep_management::BoundaryVaccuum
#define SweepIncHomoBndry chi_mesh::sweep_management::BoundaryIsotropicHomogenous
#define SweepReflectingBndry chi_mesh::sweep_management::BoundaryReflecting
#define SweepAniHeteroBndry chi_mesh::sweep_management::BoundaryIncidentHeterogeneous

#define ExceptionLocalFaceNormalsDiffer                                                            \
  std::logic_error(fname + ": Not all face normals are,"                                           \
                           " within tolerance, locally the same for the reflecting boundary"       \
                           " condition requested.")

#define ExceptionGlobalFaceNormalsDiffer                                                           \
  std::logic_error(fname + ": Not all face normals are,"                                           \
                           " within tolerance, globally the same for the reflecting boundary"      \
                           " condition requested.")

//###################################################################
/**Initializes transport related boundaries. */
void
lbs::LBSSolver::InitializeBoundaries()
{
  const std::string fname = "lbs::LBSSolver::InitializeBoundaries";
  //================================================== Determine boundary-ids
  //                                                   involved in the problem
  std::set<uint64_t> globl_unique_bids_set;
  {
    std::set<uint64_t> local_unique_bids_set;
    for (const auto& cell : grid_ptr_->local_cells)
      for (const auto& face : cell.faces_)
        if (not face.has_neighbor_) local_unique_bids_set.insert(face.neighbor_id_);

    std::vector<uint64_t> local_unique_bids(local_unique_bids_set.begin(),
                                            local_unique_bids_set.end());
    const int local_num_unique_bids = static_cast<int>(local_unique_bids.size());
    std::vector<int> recvcounts(Chi::mpi.process_count, 0);

    MPI_Allgather(&local_num_unique_bids, // sendbuf
                  1,
                  MPI_INT,           // sendcount+type
                  recvcounts.data(), // recvbuf
                  1,
                  MPI_INT,        // recvcount+type
                  Chi::mpi.comm); // comm

    std::vector<int> recvdispls(Chi::mpi.process_count, 0);

    int running_displacement = 0;
    for (int locI = 0; locI < Chi::mpi.process_count; ++locI)
    {
      recvdispls[locI] = running_displacement;
      running_displacement += recvcounts[locI];
    }

    std::vector<uint64_t> recvbuf(running_displacement, 0);

    MPI_Allgatherv(local_unique_bids.data(), // sendbuf
                   local_num_unique_bids,    // sendcount
                   MPI_UINT64_T,             // sendtype
                   recvbuf.data(),           // recvbuf
                   recvcounts.data(),        // recvcounts
                   recvdispls.data(),        // recvdispls
                   MPI_UINT64_T,             // recvtype
                   Chi::mpi.comm);           // comm

    globl_unique_bids_set = local_unique_bids_set; // give it a head start

    for (uint64_t bid : recvbuf)
      globl_unique_bids_set.insert(bid);
  }

  //================================================== Initialize default
  //                                                   incident boundary
  const size_t G = num_groups_;

  sweep_boundaries_.clear();
  for (uint64_t bid : globl_unique_bids_set)
  {
    const bool has_no_preference = boundary_preferences_.count(bid) == 0;
    const bool has_not_been_set = sweep_boundaries_.count(bid) == 0;
    if (has_no_preference and has_not_been_set)
    {
      sweep_boundaries_[bid] = mk_shrd(SweepVaccuumBndry)(G);
    } // defaulted
    else if (has_not_been_set)
    {
      const auto& bndry_pref = boundary_preferences_.at(bid);
      const auto& mg_q = bndry_pref.isotropic_mg_source;

      if (bndry_pref.type == lbs::BoundaryType::VACUUM)
        sweep_boundaries_[bid] = mk_shrd(SweepVaccuumBndry)(G);
      else if (bndry_pref.type == lbs::BoundaryType::INCIDENT_ISOTROPIC)
        sweep_boundaries_[bid] = mk_shrd(SweepIncHomoBndry)(G, mg_q);
      else if (bndry_pref.type == BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS)
      {
        sweep_boundaries_[bid] = mk_shrd(SweepAniHeteroBndry)(
          G, std::make_unique<BoundaryFunctionToLua>(bndry_pref.source_function), bid);
      }
      else if (bndry_pref.type == lbs::BoundaryType::REFLECTING)
      {
        // Locally check all faces, that subscribe to this boundary,
        // have the same normal
        typedef chi_mesh::Vector3 Vec3;
        const double EPSILON = 1.0e-12;
        std::unique_ptr<Vec3> n_ptr = nullptr;
        for (const auto& cell : grid_ptr_->local_cells)
          for (const auto& face : cell.faces_)
            if (not face.has_neighbor_ and face.neighbor_id_ == bid)
            {
              if (not n_ptr) n_ptr = std::make_unique<Vec3>(face.normal_);
              if (std::fabs(face.normal_.Dot(*n_ptr) - 1.0) > EPSILON)
                throw ExceptionLocalFaceNormalsDiffer;
            }

        // Now check globally
        const int local_has_bid = n_ptr != nullptr ? 1 : 0;
        const Vec3 local_normal = local_has_bid ? *n_ptr : Vec3(0.0, 0.0, 0.0);

        std::vector<int> locJ_has_bid(Chi::mpi.process_count, 1);
        std::vector<double> locJ_n_val(Chi::mpi.process_count * 3, 0.0);

        MPI_Allgather(&local_has_bid, // sendbuf
                      1,
                      MPI_INT,             // sendcount + datatype
                      locJ_has_bid.data(), // recvbuf
                      1,
                      MPI_INT,        // recvcount + datatype
                      Chi::mpi.comm); // communicator

        MPI_Allgather(&local_normal, // sendbuf
                      3,
                      MPI_DOUBLE,        // sendcount + datatype
                      locJ_n_val.data(), // recvbuf
                      3,
                      MPI_DOUBLE,     // recvcount + datatype
                      Chi::mpi.comm); // communicator

        Vec3 global_normal;
        for (int j = 0; j < Chi::mpi.process_count; ++j)
        {
          if (locJ_has_bid[j])
          {
            int offset = 3 * j;
            const double* n = &locJ_n_val[offset];
            const Vec3 locJ_normal(n[0], n[1], n[2]);

            if (local_has_bid)
              if (std::fabs(local_normal.Dot(locJ_normal) - 1.0) > EPSILON)
                throw ExceptionGlobalFaceNormalsDiffer;

            global_normal = locJ_normal;
          }
        }

        sweep_boundaries_[bid] = mk_shrd(SweepReflectingBndry)(
          G, global_normal, MapGeometryTypeToCoordSys(options_.geometry_type));
      }
    } // non-defaulted
  }   // for bndry id
}