#include "AAH_FLUDSCommonData.h"

#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"

#include <algorithm>

namespace chi_mesh::sweep_management
{

void
AAH_FLUDSCommonData::InitializeBetaElements(const SPDS& spds, int tag_index /*=0*/)
{
  const chi_mesh::MeshContinuum& grid = spds.Grid();
  const chi_mesh::sweep_management::SPLS& spls = spds.GetSPLS();

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // The first two major steps here are: Send delayed successor information
  // and Receive delayed predecessor information. The send portion is done
  // first because the delayed information does not follow the
  // Task Dependency Graph and hence when a location receives its delayed
  // information, the information might not have been sent yet.

  //=============================================== Send delayed successor
  // information
  const auto& location_successors = spds.GetLocationSuccessors();
  const auto& delayed_location_successors = spds.GetDelayedLocationSuccessors();
  std::vector<MPI_Request> send_requests;
  send_requests.resize(location_successors.size(), MPI_Request());
  std::vector<std::vector<int>> multi_face_indices(location_successors.size(), std::vector<int>());
  for (int deplocI = 0; deplocI < location_successors.size(); deplocI++)
  {
    int locJ = location_successors[deplocI];

    std::vector<int>::const_iterator delayed_successor =
      std::find(delayed_location_successors.begin(), delayed_location_successors.end(), locJ);
    if ((delayed_successor == delayed_location_successors.end())) continue;

    std::vector<CompactCellView> cell_views = deplocI_cell_views[deplocI];

    SerializeCellInfo(cell_views, multi_face_indices[deplocI], deplocI_face_dof_count[deplocI]);

    MPI_Isend(multi_face_indices[deplocI].data(),
              static_cast<int>(multi_face_indices[deplocI].size()),
              MPI_INT,
              locJ,
              101 + tag_index,
              Chi::mpi.comm,
              &send_requests[deplocI]);

    // TODO: Watch eager limits on sent data

    deplocI_cell_views[deplocI].clear();
    deplocI_cell_views[deplocI].shrink_to_fit();
  }

  //=============================================== Receive delayed predecessor
  //                                                information
  const auto& delayed_location_dependencies = spds.GetDelayedLocationDependencies();
  delayed_prelocI_cell_views.resize(delayed_location_dependencies.size(),
                                    std::vector<CompactCellView>());
  delayed_prelocI_face_dof_count.resize(delayed_location_dependencies.size(), 0);
  for (int prelocI = 0; prelocI < delayed_location_dependencies.size(); prelocI++)
  {
    int locJ = delayed_location_dependencies[prelocI];

    MPI_Status probe_status;
    MPI_Probe(locJ, 101 + tag_index, Chi::mpi.comm, &probe_status);

    int amount_to_receive = 0;
    MPI_Get_count(&probe_status, MPI_INT, &amount_to_receive);

    std::vector<int> face_indices;
    face_indices.resize(amount_to_receive, 0);

    MPI_Recv(face_indices.data(),
             amount_to_receive,
             MPI_INT,
             locJ,
             101 + tag_index,
             Chi::mpi.comm,
             MPI_STATUS_IGNORE);

    DeSerializeCellInfo(
      delayed_prelocI_cell_views[prelocI], &face_indices, delayed_prelocI_face_dof_count[prelocI]);
  }

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // The next two operations is to receive predecessor information followed
  // by the sending of successor information. The receives are blocking but
  // will cause a dead lock because the successor/predecessor combination
  // follows the TDG.

  //=============================================== Receive predecessor
  //                                                information
  const auto& location_dependencies = spds.GetLocationDependencies();
  prelocI_cell_views.resize(location_dependencies.size(), std::vector<CompactCellView>());
  prelocI_face_dof_count.resize(location_dependencies.size(), 0);
  for (int prelocI = 0; prelocI < location_dependencies.size(); prelocI++)
  {
    int locJ = location_dependencies[prelocI];

    MPI_Status probe_status;
    MPI_Probe(locJ, 101 + tag_index, Chi::mpi.comm, &probe_status);

    int amount_to_receive = 0;
    MPI_Get_count(&probe_status, MPI_INT, &amount_to_receive);

    std::vector<int> face_indices;
    face_indices.resize(amount_to_receive, 0);

    MPI_Recv(face_indices.data(),
             amount_to_receive,
             MPI_INT,
             locJ,
             101 + tag_index,
             Chi::mpi.comm,
             MPI_STATUS_IGNORE);

    DeSerializeCellInfo(
      prelocI_cell_views[prelocI], &face_indices, prelocI_face_dof_count[prelocI]);
  }

  //=============================================== Send successor information
  for (int deplocI = 0; deplocI < location_successors.size(); deplocI++)
  {
    int locJ = location_successors[deplocI];

    auto delayed_successor =
      std::find(delayed_location_successors.begin(), delayed_location_successors.end(), locJ);
    if ((delayed_successor != delayed_location_successors.end())) continue;

    std::vector<CompactCellView> cell_views = deplocI_cell_views[deplocI];

    SerializeCellInfo(cell_views, multi_face_indices[deplocI], deplocI_face_dof_count[deplocI]);

    MPI_Isend(multi_face_indices[deplocI].data(),
              static_cast<int>(multi_face_indices[deplocI].size()),
              MPI_INT,
              locJ,
              101 + tag_index,
              Chi::mpi.comm,
              &send_requests[deplocI]);

    // TODO: Watch eager limits on sent data

    deplocI_cell_views[deplocI].clear();
    deplocI_cell_views[deplocI].shrink_to_fit();
  }

  //================================================== Verify sends completed
  for (int deplocI = 0; deplocI < location_successors.size(); deplocI++)
    MPI_Wait(&send_requests[deplocI], MPI_STATUS_IGNORE);
  multi_face_indices.clear();
  multi_face_indices.shrink_to_fit();

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  // In the next process we loop over cells in the sweep order and perform
  // the non-local face mappings. This is dependent on having the compact
  // cellviews on the partition interfaces.

  //================================================== Loop over cells in sorder
  for (int csoi = 0; csoi < spls.item_id.size(); csoi++)
  {
    int cell_local_index = spls.item_id[csoi];
    const auto& cell = grid.local_cells[cell_local_index];

    NonLocalIncidentMapping(cell, spds);
  } // for csoi

  deplocI_cell_views.clear();
  deplocI_cell_views.shrink_to_fit();

  prelocI_cell_views.clear();
  prelocI_cell_views.shrink_to_fit();

  delayed_prelocI_cell_views.clear();
  delayed_prelocI_cell_views.shrink_to_fit();

  //================================================== Clear unneccesary data
  auto empty_vector = std::vector<std::vector<CompactCellView>>(0);
  deplocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  prelocI_cell_views.swap(empty_vector);

  empty_vector = std::vector<std::vector<CompactCellView>>(0);
  delayed_prelocI_cell_views.swap(empty_vector);
}

} // namespace chi_mesh::sweep_management