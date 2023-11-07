#include "framework/mesh/SweepUtilities/FLUDS/AAH_FLUDS.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/math/chi_math.h"

namespace chi_mesh::sweep_management
{

/**This constructor initializes an auxiliary FLUDS based
 * on a primary FLUDS. The restriction here is that the
 * auxiliary FLUDS has the exact same sweep ordering as the
 * primary FLUDS.*/
AAH_FLUDS::AAH_FLUDS(size_t num_groups, size_t num_angles, const AAH_FLUDSCommonData& common_data)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()), common_data_(common_data)
{
  // Adjusting for different group aggregate
  for (auto& val : common_data_.local_psi_n_block_stride)
    local_psi_Gn_block_strideG.push_back(val * num_groups_);

  delayed_local_psi_Gn_block_strideG = common_data_.delayed_local_psi_Gn_block_stride * num_groups_;
}

/**Given a sweep ordering index, the outgoing face counter,
 * the outgoing face dof, this function computes the location
 * of this position's upwind psi in the local upwind psi vector
 * and returns a reference to it.*/
double*
AAH_FLUDS::OutgoingPsi(int cell_so_index, int outb_face_counter, int face_dof, int n)
{
  // Face category
  int fc = common_data_.so_cell_outb_face_face_category[cell_so_index][outb_face_counter];

  if (fc >= 0)
  {
    size_t index = local_psi_Gn_block_strideG[fc] * n +
                   common_data_.so_cell_outb_face_slot_indices[cell_so_index][outb_face_counter] *
                     common_data_.local_psi_stride[fc] * num_groups_ +
                   face_dof * num_groups_;

    return &local_psi_[fc][index];
  }
  else
  {
    size_t index = delayed_local_psi_Gn_block_strideG * n +
                   common_data_.so_cell_outb_face_slot_indices[cell_so_index][outb_face_counter] *
                     common_data_.delayed_local_psi_stride * num_groups_ +
                   face_dof * num_groups_;

    return &delayed_local_psi_[index];
  }
}

/**Given a outbound face counter this method returns a pointer
 * to the location*/
double*
AAH_FLUDS::NLOutgoingPsi(int outb_face_counter, int face_dof, int n)
{
  if (outb_face_counter > common_data_.nonlocal_outb_face_deplocI_slot.size())
  {
    Chi::log.LogAllError() << "Invalid number of outb_face_counter " << outb_face_counter
                           << " max allowed "
                           << common_data_.nonlocal_outb_face_deplocI_slot.size();
    Chi::Exit(EXIT_FAILURE);
  }

  int depLocI = common_data_.nonlocal_outb_face_deplocI_slot[outb_face_counter].first;
  int slot = common_data_.nonlocal_outb_face_deplocI_slot[outb_face_counter].second;
  int nonlocal_psi_Gn_blockstride = common_data_.deplocI_face_dof_count[depLocI];

  int index =
    nonlocal_psi_Gn_blockstride * num_groups_ * n + slot * num_groups_ + face_dof * num_groups_;

  if ((index < 0) || (index > deplocI_outgoing_psi_[depLocI].size()))
  {
    Chi::log.LogAllError() << "Invalid index " << index << " encountered in non-local outgoing Psi"
                           << " max allowed " << deplocI_outgoing_psi_[depLocI].size();
    Chi::Exit(EXIT_FAILURE);
  }

  return &deplocI_outgoing_psi_[depLocI][index];
}

/**Given a sweep ordering index, the incoming face counter,
 * the incoming face dof, this function computes the location
 * where to store this position's outgoing psi and returns a reference
 * to it.*/
double*
AAH_FLUDS::UpwindPsi(int cell_so_index, int inc_face_counter, int face_dof, int g, int n)
{
  // Face category
  int fc = common_data_.so_cell_inco_face_face_category[cell_so_index][inc_face_counter];

  if (fc >= 0)
  {
    size_t index =
      local_psi_Gn_block_strideG[fc] * n +
      common_data_.so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].slot_address *
        common_data_.local_psi_stride[fc] * num_groups_ +
      common_data_.so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter]
          .upwind_dof_mapping[face_dof] *
        num_groups_ +
      g;

    return &local_psi_[fc][index];
  }
  else
  {
    size_t index =
      delayed_local_psi_Gn_block_strideG * n +
      common_data_.so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter].slot_address *
        common_data_.delayed_local_psi_stride * num_groups_ +
      common_data_.so_cell_inco_face_dof_indices[cell_so_index][inc_face_counter]
          .upwind_dof_mapping[face_dof] *
        num_groups_ +
      g;

    return &delayed_local_psi_old_[index];
  }
}

/**Given a sweep ordering index, the incoming face counter,
 * the incoming face dof, this function computes the location
 * where to obtain the position's upwind psi.*/
double*
AAH_FLUDS::NLUpwindPsi(int nonl_inc_face_counter, int face_dof, int g, int n)
{
  int prelocI = common_data_.nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].first;

  if (prelocI >= 0)
  {
    int nonlocal_psi_Gn_blockstride = common_data_.prelocI_face_dof_count[prelocI];
    int slot = common_data_.nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].second.first;

    int mapped_dof = common_data_.nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter]
                       .second.second[face_dof];

    int index = nonlocal_psi_Gn_blockstride * num_groups_ * n + slot * num_groups_ +
                mapped_dof * num_groups_ + g;

    return &prelocI_outgoing_psi_[prelocI][index];
  }
  else
  {
    prelocI = common_data_.delayed_nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].first;

    int nonlocal_psi_Gn_blockstride = common_data_.delayed_prelocI_face_dof_count[prelocI];
    int slot =
      common_data_.delayed_nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter].second.first;

    int mapped_dof = common_data_.delayed_nonlocal_inc_face_prelocI_slot_dof[nonl_inc_face_counter]
                       .second.second[face_dof];

    int index = nonlocal_psi_Gn_blockstride * num_groups_ * n + slot * num_groups_ +
                mapped_dof * num_groups_ + g;

    return &delayed_prelocI_outgoing_psi_old_[prelocI][index];
  }
}

size_t
AAH_FLUDS::GetPrelocIFaceDOFCount(int prelocI) const
{
  return common_data_.prelocI_face_dof_count[prelocI];
}
size_t
AAH_FLUDS::GetDelayedPrelocIFaceDOFCount(int prelocI) const
{
  return common_data_.delayed_prelocI_face_dof_count[prelocI];
}
size_t
AAH_FLUDS::GetDeplocIFaceDOFCount(int deplocI) const
{
  return common_data_.deplocI_face_dof_count[deplocI];
}

void
AAH_FLUDS::ClearLocalAndReceivePsi()
{
  auto empty_vector = std::vector<std::vector<double>>(0);
  local_psi_.swap(empty_vector);

  empty_vector = std::vector<std::vector<double>>(0);
  prelocI_outgoing_psi_.swap(empty_vector);
}

void
AAH_FLUDS::ClearSendPsi()
{
  deplocI_outgoing_psi_.clear();
}

void
AAH_FLUDS::AllocateInternalLocalPsi(size_t num_grps, size_t num_angles)
{
  local_psi_.resize(common_data_.num_face_categories);
  // fc = face category
  for (size_t fc = 0; fc < common_data_.num_face_categories; fc++)
  {
    local_psi_[fc].resize(common_data_.local_psi_stride[fc] *
                            common_data_.local_psi_max_elements[fc] * num_grps * num_angles,
                          0.0);
  }
}

void
AAH_FLUDS::AllocateOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_sucs)
{
  deplocI_outgoing_psi_.resize(num_loc_sucs, std::vector<double>());
  for (size_t deplocI = 0; deplocI < num_loc_sucs; deplocI++)
  {
    deplocI_outgoing_psi_[deplocI].resize(
      common_data_.deplocI_face_dof_count[deplocI] * num_grps * num_angles, 0.0);
  }
}

void
AAH_FLUDS::AllocateDelayedLocalPsi(size_t num_grps, size_t num_angles)
{
  delayed_local_psi_.resize(common_data_.delayed_local_psi_stride *
                              common_data_.delayed_local_psi_max_elements * num_grps * num_angles,
                            0.0);

  delayed_local_psi_old_.resize(common_data_.delayed_local_psi_stride *
                                  common_data_.delayed_local_psi_max_elements * num_grps *
                                  num_angles,
                                0.0);
}

void
AAH_FLUDS::AllocatePrelocIOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_deps)
{
  prelocI_outgoing_psi_.resize(num_loc_deps, std::vector<double>());
  for (size_t prelocI = 0; prelocI < num_loc_deps; prelocI++)
  {
    prelocI_outgoing_psi_[prelocI].resize(
      common_data_.prelocI_face_dof_count[prelocI] * num_grps * num_angles, 0.0);
  }
}

void
AAH_FLUDS::AllocateDelayedPrelocIOutgoingPsi(size_t num_grps,
                                             size_t num_angles,
                                             size_t num_loc_deps)
{
  delayed_prelocI_outgoing_psi_.clear();
  delayed_prelocI_outgoing_psi_.resize(num_loc_deps);

  delayed_prelocI_outgoing_psi_old_.clear();
  delayed_prelocI_outgoing_psi_old_.resize(num_loc_deps);

  for (int prelocI = 0; prelocI < num_loc_deps; prelocI++)
  {
    const int num_nodes = common_data_.delayed_prelocI_face_dof_count[prelocI];

    uint64_t buff_size = num_nodes * num_grps * num_angles;

    delayed_prelocI_outgoing_psi_[prelocI].resize(buff_size, 0.0);
    delayed_prelocI_outgoing_psi_old_[prelocI].resize(buff_size, 0.0);
  }
}

std::vector<double>&
AAH_FLUDS::DelayedLocalPsi()
{
  return delayed_local_psi_;
}

std::vector<double>&
AAH_FLUDS::DelayedLocalPsiOld()
{
  return delayed_local_psi_old_;
}

std::vector<std::vector<double>>&
AAH_FLUDS::DeplocIOutgoingPsi()
{
  return deplocI_outgoing_psi_;
}

std::vector<std::vector<double>>&
AAH_FLUDS::PrelocIOutgoingPsi()
{
  return prelocI_outgoing_psi_;
}

std::vector<std::vector<double>>&
AAH_FLUDS::DelayedPrelocIOutgoingPsi()
{
  return delayed_prelocI_outgoing_psi_;
}
std::vector<std::vector<double>>&
AAH_FLUDS::DelayedPrelocIOutgoingPsiOld()
{
  return delayed_prelocI_outgoing_psi_old_;
}

} // namespace chi_mesh::sweep_management
