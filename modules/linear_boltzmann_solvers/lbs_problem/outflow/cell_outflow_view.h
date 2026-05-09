// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace opensn
{

/**
 * Non-owning view of group-wise outflow values for the faces of one local cell.
 *
 * Face offsets are relative to the first outflow value assigned to the cell. A negative face offset
 * marks a face that does not contribute outflow storage.
 */
class CellOutflowView
{
public:
  CellOutflowView() = default;

  /**
   * Construct a view for a cell with a fixed face count and energy group count.
   * \param num_faces Number of faces on the viewed cell.
   * \param num_groups Number of energy groups stored per boundary face.
   */
  CellOutflowView(std::size_t num_faces, unsigned int num_groups)
    : num_faces_(num_faces), num_groups_(num_groups)
  {
  }

  /// Return whether the view was initialized with at least one face or group.
  bool IsInitialized() const { return num_faces_ > 0 or num_groups_ > 0; }

  /// Return whether boundary-face offset metadata has been allocated.
  bool HasFaceOffsets() const { return not face_offsets_.empty(); }

  /// Allocate boundary-face offset metadata.
  void InitializeFaceOffsets()
  {
    face_offsets_.reserve(num_faces_);
    face_offsets_.assign(num_faces_, static_cast<std::int64_t>(-1));
  }

  /**
   * Assign the first outflow value for this cell.
   * \param outflow Pointer to the first group-wise outflow value owned by the outflow bank.
   * \note The view does not own the assigned storage.
   */
  void Assign(double* outflow) { outflow_ = outflow; }

  /**
   * Set the group-wise outflow offset for a face.
   * \param f Face index local to the viewed cell.
   * \param offset Offset of the face's first group value relative to the cell outflow pointer.
   */
  void SetFaceOffset(std::size_t f, std::int64_t offset)
  {
    assert(f < num_faces_ && "CellOutflowView::SetFaceOffset face index out of range.");
    assert(f < face_offsets_.size() &&
           "CellOutflowView::SetFaceOffset called before face offsets were initialized.");
    face_offsets_[f] = offset;
  }

  /**
   * Set one face-group outflow value to zero when storage exists.
   * \param f Face index local to the viewed cell.
   * \param g Energy group index.
   */
  void Zero(std::size_t f, unsigned int g)
  {
    const auto offset = GetOffset(f, g);
    if (offset >= 0)
      outflow_[offset] = 0.0;
  }

  /**
   * Add a contribution to one face-group outflow value when storage exists.
   * \param f Face index local to the viewed cell.
   * \param g Energy group index.
   * \param intS_mu_psi Surface-integral contribution to the outflow tally.
   */
  void Add(std::size_t f, unsigned int g, double intS_mu_psi)
  {
    const auto offset = GetOffset(f, g);
    if (offset >= 0)
      outflow_[offset] += intS_mu_psi;
  }

  /**
   * Return one face-group outflow value, or zero when storage is absent.
   * \param f Face index local to the viewed cell.
   * \param g Energy group index.
   * \return Stored outflow tally for the face and group.
   */
  double Get(std::size_t f, unsigned int g) const
  {
    const auto offset = GetOffset(f, g);
    return offset >= 0 ? outflow_[offset] : 0.0;
  }

private:
  /**
   * Return the storage offset for a face-group pair, or a negative value when storage is absent.
   * \param f Face index local to the viewed cell.
   * \param g Energy group index.
   * \return Relative offset into the assigned outflow storage.
   */
  std::int64_t GetOffset(std::size_t f, unsigned int g) const
  {
    if (f >= num_faces_)
      return -1;

    assert(g < num_groups_ && "CellOutflowView group index out of range.");
    if (outflow_ == nullptr or f >= face_offsets_.size() or face_offsets_[f] < 0 or
        g >= num_groups_)
      return -1;

    return face_offsets_[f] + static_cast<std::int64_t>(g);
  }

  /// Number of faces on the viewed cell.
  std::size_t num_faces_ = 0;

  /// Number of energy groups stored per boundary face.
  unsigned int num_groups_ = 0;

  /// Pointer to the first outflow value for the viewed cell.
  double* outflow_ = nullptr;

  /**
   * Relative offset into the cell outflow array for each face.
   * Non-boundary faces are assigned to -1.
   */
  std::vector<std::int64_t> face_offsets_;
};

} // namespace opensn
