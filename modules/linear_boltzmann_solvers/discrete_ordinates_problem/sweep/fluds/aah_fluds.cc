// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/logging/log.h"
#include "framework/math/math.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

static inline void
UpdateRange(std::vector<std::vector<double>>& target, std::vector<std::span<double>>& view)
{
  view.resize(target.size());
  for (std::size_t i = 0; i < target.size(); ++i)
  {
    view[i] = std::span<double>(target[i]);
  }
}

AAH_FLUDS::AAH_FLUDS(unsigned int num_groups,
                     size_t num_angles,
                     const AAH_FLUDSCommonData& common_data)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    delayed_local_psi_Gn_block_strideG_(common_data_.delayed_local_psi_Gn_block_stride_ *
                                        num_groups_)
{
  CALI_CXX_MARK_SCOPE("AAH_FLUDS::AAH_FLUDS");

  // Adjusting for different group aggregate
  for (const auto& val : common_data_.local_psi_n_block_stride_)
    local_psi_Gn_block_strideG_.push_back(val * num_groups_);
}

double*
AAH_FLUDS::OutgoingPsi(std::size_t cell_so_index,
                       int outb_face_counter,
                       std::size_t face_dof,
                       std::size_t n)
{
  // Face category
  int fc = common_data_.so_cell_outb_face_face_category_[cell_so_index][outb_face_counter];

  if (fc >= 0)
  {
    size_t index = local_psi_Gn_block_strideG_[fc] * n +
                   common_data_.so_cell_outb_face_slot_indices_[cell_so_index][outb_face_counter] *
                     common_data_.local_psi_stride_[fc] * num_groups_ +
                   face_dof * num_groups_;

    return &local_psi_[fc][index];
  }
  else
  {
    size_t index = delayed_local_psi_Gn_block_strideG_ * n +
                   common_data_.so_cell_outb_face_slot_indices_[cell_so_index][outb_face_counter] *
                     common_data_.delayed_local_psi_stride_ * num_groups_ +
                   face_dof * num_groups_;

    return &delayed_local_psi_[index];
  }
}

double*
AAH_FLUDS::NLOutgoingPsi(int outb_face_counter, std::size_t face_dof, std::size_t n)
{
  if (outb_face_counter > common_data_.nonlocal_outb_face_deplocI_slot_.size())
  {
    std::ostringstream oss;
    oss << "AAH_FLUDS: Invalid value for outb_face_counter " << outb_face_counter << " (max "
        << "allowed = " << common_data_.nonlocal_outb_face_deplocI_slot_.size() << ")";
    throw std::runtime_error(oss.str());
  }

  int depLocI = common_data_.nonlocal_outb_face_deplocI_slot_[outb_face_counter].first;
  int slot = common_data_.nonlocal_outb_face_deplocI_slot_[outb_face_counter].second;
  auto nonlocal_psi_Gn_blockstride = common_data_.deplocI_face_dof_count_[depLocI];

  std::size_t index = nonlocal_psi_Gn_blockstride * num_groups_ * n +
                      static_cast<std::size_t>(slot) * num_groups_ + face_dof * num_groups_;

  if ((index < 0) or (index > deplocI_outgoing_psi_[depLocI].size()))
  {
    std::stringstream oss;
    oss << "AAH_FLUDS: Invalid index " << index
        << " encountered in non-local outgoing psi (max allowed = "
        << deplocI_outgoing_psi_[depLocI].size() << ")";
    throw std::runtime_error(oss.str());
  }

  return &deplocI_outgoing_psi_[depLocI][index];
}

double*
AAH_FLUDS::UpwindPsi(std::size_t cell_so_index,
                     int inc_face_counter,
                     std::size_t face_dof,
                     unsigned int g,
                     std::size_t n)
{
  // Face category
  int fc = common_data_.so_cell_inco_face_face_category_[cell_so_index][inc_face_counter];

  if (fc >= 0)
  {
    size_t index =
      local_psi_Gn_block_strideG_[fc] * n +
      common_data_.so_cell_inco_face_dof_indices_[cell_so_index][inc_face_counter].slot_address *
        common_data_.local_psi_stride_[fc] * num_groups_ +
      static_cast<size_t>(
        common_data_.so_cell_inco_face_dof_indices_[cell_so_index][inc_face_counter]
          .upwind_dof_mapping[face_dof]) *
        num_groups_ +
      g;

    return &local_psi_[fc][index];
  }
  else
  {
    size_t index =
      delayed_local_psi_Gn_block_strideG_ * n +
      common_data_.so_cell_inco_face_dof_indices_[cell_so_index][inc_face_counter].slot_address *
        common_data_.delayed_local_psi_stride_ * num_groups_ +
      static_cast<size_t>(
        common_data_.so_cell_inco_face_dof_indices_[cell_so_index][inc_face_counter]
          .upwind_dof_mapping[face_dof]) *
        num_groups_ +
      g;

    return &delayed_local_psi_old_[index];
  }
}

double*
AAH_FLUDS::NLUpwindPsi(int nonl_inc_face_counter,
                       std::size_t face_dof,
                       unsigned int g,
                       std::size_t n)
{
  int prelocI = common_data_.nonlocal_inc_face_prelocI_slot_dof_[nonl_inc_face_counter].first;

  if (prelocI >= 0)
  {
    std::size_t nonlocal_psi_Gn_blockstride = common_data_.prelocI_face_dof_count_[prelocI];
    auto slot =
      common_data_.nonlocal_inc_face_prelocI_slot_dof_[nonl_inc_face_counter].second.first;

    auto mapped_dof = common_data_.nonlocal_inc_face_prelocI_slot_dof_[nonl_inc_face_counter]
                        .second.second[face_dof];

    std::size_t index = nonlocal_psi_Gn_blockstride * num_groups_ * n + slot * num_groups_ +
                        mapped_dof * num_groups_ + g;

    return &prelocI_outgoing_psi_[prelocI][index];
  }
  else
  {
    prelocI = common_data_.delayed_nonlocal_inc_face_prelocI_slot_dof_[nonl_inc_face_counter].first;

    std::size_t nonlocal_psi_Gn_blockstride = common_data_.delayed_prelocI_face_dof_count_[prelocI];
    auto slot =
      common_data_.delayed_nonlocal_inc_face_prelocI_slot_dof_[nonl_inc_face_counter].second.first;

    auto mapped_dof =
      common_data_.delayed_nonlocal_inc_face_prelocI_slot_dof_[nonl_inc_face_counter]
        .second.second[face_dof];

    std::size_t index = nonlocal_psi_Gn_blockstride * num_groups_ * n + slot * num_groups_ +
                        mapped_dof * num_groups_ + g;

    return &delayed_prelocI_outgoing_psi_old_[prelocI][index];
  }
}

size_t
AAH_FLUDS::GetPrelocIFaceDOFCount(std::size_t prelocI) const
{
  return common_data_.prelocI_face_dof_count_[prelocI];
}

size_t
AAH_FLUDS::GetDelayedPrelocIFaceDOFCount(std::size_t prelocI) const
{
  return common_data_.delayed_prelocI_face_dof_count_[prelocI];
}

size_t
AAH_FLUDS::GetDeplocIFaceDOFCount(std::size_t deplocI) const
{
  return common_data_.deplocI_face_dof_count_[deplocI];
}

void
AAH_FLUDS::ClearLocalAndReceivePsi()
{
  local_psi_.clear();
  prelocI_outgoing_psi_.clear();
  prelocI_outgoing_psi_view_.clear();
}

void
AAH_FLUDS::ClearSendPsi()
{
  deplocI_outgoing_psi_.clear();
  deplocI_outgoing_psi_view_.clear();
}

void
AAH_FLUDS::AllocateInternalLocalPsi()
{
  local_psi_.resize(common_data_.num_face_categories_);
  // fc = face category
  for (size_t fc = 0; fc < common_data_.num_face_categories_; ++fc)
  {
    local_psi_[fc].resize(common_data_.local_psi_stride_[fc] *
                            common_data_.local_psi_max_elements_[fc] * num_groups_and_angles_,
                          0.0);
  }
}

void
AAH_FLUDS::AllocateOutgoingPsi()
{
  std::size_t num_loc_sucs = spds_.GetLocationSuccessors().size();
  deplocI_outgoing_psi_.resize(num_loc_sucs, std::vector<double>());
  for (size_t deplocI = 0; deplocI < num_loc_sucs; ++deplocI)
  {
    deplocI_outgoing_psi_[deplocI].resize(
      common_data_.deplocI_face_dof_count_[deplocI] * num_groups_and_angles_, 0.0);
  }
  UpdateRange(deplocI_outgoing_psi_, deplocI_outgoing_psi_view_);
}

void
AAH_FLUDS::AllocateDelayedLocalPsi()
{
  std::size_t delayed_local_psi_size = common_data_.delayed_local_psi_stride_ *
                                       common_data_.delayed_local_psi_max_elements_ *
                                       num_groups_and_angles_;
  delayed_local_psi_.resize(delayed_local_psi_size, 0.0);
  delayed_local_psi_view_ = std::span<double>(delayed_local_psi_);

  delayed_local_psi_old_.resize(delayed_local_psi_size, 0.0);
  delayed_local_psi_old_view_ = std::span<double>(delayed_local_psi_old_);
}

void
AAH_FLUDS::AllocatePrelocIOutgoingPsi()
{
  std::size_t num_loc_deps = spds_.GetLocationDependencies().size();
  prelocI_outgoing_psi_.resize(num_loc_deps, std::vector<double>());
  for (size_t prelocI = 0; prelocI < num_loc_deps; ++prelocI)
  {
    prelocI_outgoing_psi_[prelocI].resize(
      common_data_.prelocI_face_dof_count_[prelocI] * num_groups_and_angles_, 0.0);
  }
  UpdateRange(prelocI_outgoing_psi_, prelocI_outgoing_psi_view_);
}

void
AAH_FLUDS::AllocateDelayedPrelocIOutgoingPsi()
{
  std::size_t num_loc_deps = spds_.GetDelayedLocationDependencies().size();
  delayed_prelocI_outgoing_psi_.resize(num_loc_deps);
  delayed_prelocI_outgoing_psi_old_.resize(num_loc_deps);

  for (size_t prelocI = 0; prelocI < num_loc_deps; ++prelocI)
  {
    const auto num_nodes = common_data_.delayed_prelocI_face_dof_count_[prelocI];

    uint64_t buff_size = num_nodes * num_groups_and_angles_;

    delayed_prelocI_outgoing_psi_[prelocI].resize(buff_size, 0.0);
    delayed_prelocI_outgoing_psi_old_[prelocI].resize(buff_size, 0.0);
  }
  UpdateRange(delayed_prelocI_outgoing_psi_, delayed_prelocI_outgoing_psi_view_);
  UpdateRange(delayed_prelocI_outgoing_psi_old_, delayed_prelocI_outgoing_psi_old_view_);
}

void
AAH_FLUDS::SetDelayedOutgoingPsiNewToOld()
{
  delayed_prelocI_outgoing_psi_old_ = delayed_prelocI_outgoing_psi_;
  UpdateRange(delayed_prelocI_outgoing_psi_old_, delayed_prelocI_outgoing_psi_old_view_);
}

void
AAH_FLUDS::SetDelayedOutgoingPsiOldToNew()
{
  delayed_prelocI_outgoing_psi_ = delayed_prelocI_outgoing_psi_old_;
  UpdateRange(delayed_prelocI_outgoing_psi_, delayed_prelocI_outgoing_psi_view_);
}

void
AAH_FLUDS::SetDelayedLocalPsiNewToOld()
{
  delayed_local_psi_old_ = delayed_local_psi_;
  delayed_local_psi_old_view_ = std::span<double>(delayed_local_psi_old_);
}

void
AAH_FLUDS::SetDelayedLocalPsiOldToNew()
{
  delayed_local_psi_ = delayed_local_psi_old_;
  delayed_local_psi_view_ = std::span<double>(delayed_local_psi_);
}

} // namespace opensn
