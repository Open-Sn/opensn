// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#ifdef __HIPCC__
#include <hip/hip_version.h>
#endif
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace opensn::gpu_kernel
{

/// Bytes reserved for non-metadata arguments and alignment, checked in cbcd_fused_main.h.
inline constexpr std::size_t CBCD_OTHER_PARAMETER_ALLOWANCE = 64;

inline constexpr std::size_t CBCD_SMALL_PARAMETER_BYTES = 4096;
inline constexpr std::size_t CBCD_HIP_PARAMETER_BYTES = 32768;
inline constexpr std::size_t CBCD_SMALL_LAUNCH_BYTES =
  CBCD_SMALL_PARAMETER_BYTES - CBCD_OTHER_PARAMETER_ALLOWANCE;

inline constexpr std::size_t CBCD_MEDIUM_LAUNCH_BYTES = 16384;

/// Launch metadata byte budgets by compiler profile. CUDA uses the larger parameter budget
/// where supported. See the CUDA C++ Programming Guide, “__global__ Function Parameters.”
/// HIP uses a fixed budget. SYCL passes a pointer to device USM.
inline constexpr std::size_t CBCD_LAUNCH_BYTES =
#ifdef __NVCC__
  (__CUDACC_VER_MAJOR__ > 12 or (__CUDACC_VER_MAJOR__ == 12 and __CUDACC_VER_MINOR__ >= 1))
    ? CBCD_MEDIUM_LAUNCH_BYTES
    : CBCD_SMALL_LAUNCH_BYTES;
#elif defined(__HIPCC__)
  (HIP_VERSION_MAJOR > 7 or (HIP_VERSION_MAJOR == 7 and HIP_VERSION_MINOR >= 2))
    ? CBCD_HIP_PARAMETER_BYTES - CBCD_OTHER_PARAMETER_ALLOWANCE
    : CBCD_SMALL_LAUNCH_BYTES;
#else
  CBCD_SMALL_LAUNCH_BYTES;
#endif

/// Metadata size classes in bytes, smallest first. Dispatch uses the smallest distinct class that
/// fits because the whole declared extent of a by-value parameter is copied.
inline constexpr std::array<std::size_t, 3> CBCD_LAUNCH_CLASS_BYTES{
  std::min(CBCD_SMALL_LAUNCH_BYTES, CBCD_LAUNCH_BYTES),
  std::min(CBCD_MEDIUM_LAUNCH_BYTES, CBCD_LAUNCH_BYTES),
  CBCD_LAUNCH_BYTES};

/// Number of positions in the launch-function table, including duplicate size classes.
inline constexpr std::size_t CBCD_NUM_LAUNCH_CLASSES = CBCD_LAUNCH_CLASS_BYTES.size();

/// Choose how much of the metadata extent to reserve for angle-set slots.
constexpr std::uint32_t
CBCDMaxBatchSlots(std::size_t bytes)
{
  return bytes >= CBCD_MEDIUM_LAUNCH_BYTES ? 64 : (bytes >= CBCD_SMALL_LAUNCH_BYTES ? 32 : 16);
}

/// Number of uint32 counts before a launch's slots: used slots and used entries.
inline constexpr std::size_t CBCD_LAUNCH_HEADER_ENTRIES = 2;

inline constexpr auto CBCD_CELL_ID_BITS = std::numeric_limits<std::uint32_t>::digits;
inline constexpr auto CBCD_FACE_OFFSET_BITS = std::numeric_limits<std::uint16_t>::digits;
inline constexpr auto CBCD_OUTGOING_OFFSET_SHIFT = CBCD_CELL_ID_BITS + CBCD_FACE_OFFSET_BITS;
inline constexpr std::uint64_t CBCD_FACE_OFFSET_MASK = std::numeric_limits<std::uint16_t>::max();

/// One angle set's ready cells in a fused CBCD launch.
struct CBCD_BatchSlot
{
  /// One past the last thread block of this slot within the launch.
  std::uint32_t block_end;
  /// Index in CBCD_Launch::entries of the slot's first cell; its payload offsets follow its cells.
  std::uint32_t first_entry;
  std::uint32_t num_cells;
  /// Number of payload offsets of incoming nonlocal faces; those of outgoing ones follow.
  std::uint32_t num_incoming_faces;
  std::uint32_t angle_set_id;
  std::uint16_t block_size_x;
  std::uint16_t num_stride_blocks;
};

/// Fused launch metadata bounded by bytes. Each slot stores cell entries, then incoming and
/// outgoing face offsets in doubles, in each cell's face order. A cell entry packs its local ID
/// (bits 0–31) and its first incoming (bits 32–47) and outgoing (bits 48–63) offset indices within
/// the slot.
template <std::size_t bytes>
struct CBCD_LaunchOf
{
  static constexpr std::uint32_t max_slots = CBCDMaxBatchSlots(bytes);
  static_assert(bytes >= CBCD_LAUNCH_HEADER_ENTRIES * sizeof(std::uint32_t) +
                           max_slots * sizeof(CBCD_BatchSlot),
                "Device CBC launch extent is too small for its header and batch slots.");
  static constexpr std::uint32_t max_entries =
    static_cast<std::uint32_t>((bytes - CBCD_LAUNCH_HEADER_ENTRIES * sizeof(std::uint32_t) -
                                max_slots * sizeof(CBCD_BatchSlot)) /
                               sizeof(std::uint64_t));

  static constexpr std::uint64_t EncodeCellData(std::uint32_t cell_local_id,
                                                std::uint32_t first_incoming,
                                                std::uint32_t first_outgoing)
  {
    return std::uint64_t{cell_local_id} | (std::uint64_t{first_incoming} << CBCD_CELL_ID_BITS) |
           (std::uint64_t{first_outgoing} << CBCD_OUTGOING_OFFSET_SHIFT);
  }

  static constexpr bool Fits(std::uint32_t num_slots, std::uint32_t num_entries)
  {
    return num_slots <= max_slots and num_entries <= max_entries;
  }

  /// Copy used metadata. The source counts must fit the destination's capacities.
  template <std::size_t other_bytes>
  void CopyUsed(const CBCD_LaunchOf<other_bytes>& launch)
  {
    num_slots = launch.num_slots;
    num_entries = launch.num_entries;
    std::copy_n(launch.slots.begin(), num_slots, slots.begin());
    std::copy_n(launch.entries.begin(), num_entries, entries.begin());
  }

  std::uint32_t num_slots;
  std::uint32_t num_entries;
  /// Slots, ordered by block range.
  std::array<CBCD_BatchSlot, max_slots> slots;
  std::array<std::uint64_t, max_entries> entries;
};

template <std::size_t launch_class>
using CBCD_LaunchClass = CBCD_LaunchOf<CBCD_LAUNCH_CLASS_BYTES[launch_class]>;

// Launch of the full parameter budget, in which the host builds every launch.
using CBCD_Launch = CBCD_LaunchClass<CBCD_NUM_LAUNCH_CLASSES - 1>;

/// Skip classes duplicated by the next class. The full-budget class is always used.
constexpr bool
IsCBCDLaunchClassUsed(std::size_t launch_class)
{
  return launch_class + 1 == CBCD_NUM_LAUNCH_CLASSES or
         CBCD_LAUNCH_CLASS_BYTES[launch_class] < CBCD_LAUNCH_CLASS_BYTES[launch_class + 1];
}

/// Find the smallest used class that fits both counts. Both must fit the full class.
template <std::size_t launch_class = 0>
constexpr std::size_t
GetCBCDLaunchClass(std::uint32_t num_slots, std::uint32_t num_entries)
{
  if constexpr (launch_class + 1 == CBCD_NUM_LAUNCH_CLASSES)
    return launch_class; // NOLINT(bugprone-branch-clone)
  else if (IsCBCDLaunchClassUsed(launch_class) and
           CBCD_LaunchClass<launch_class>::Fits(num_slots, num_entries))
    return launch_class;
  else
    return GetCBCDLaunchClass<launch_class + 1>(num_slots, num_entries);
}

template <std::size_t bytes>
constexpr bool
IsCBCDLaunchLayoutValid()
{
  using Launch = CBCD_LaunchOf<bytes>;
  return std::is_standard_layout_v<Launch> and std::is_trivially_copyable_v<Launch> and
         offsetof(Launch, num_slots) == 0 and
         offsetof(Launch, num_entries) == sizeof(std::uint32_t) and
         offsetof(Launch, slots) == CBCD_LAUNCH_HEADER_ENTRIES * sizeof(std::uint32_t) and
         offsetof(Launch, entries) == CBCD_LAUNCH_HEADER_ENTRIES * sizeof(std::uint32_t) +
                                        Launch::max_slots * sizeof(CBCD_BatchSlot) and
         sizeof(Launch) == bytes and alignof(Launch) == alignof(std::uint64_t) and
         Launch::max_entries <= CBCD_FACE_OFFSET_MASK;
}

// The wire/parameter layout is five uint32 fields followed by two uint16 fields.
static_assert(std::is_standard_layout_v<CBCD_BatchSlot> and
                std::is_trivially_copyable_v<CBCD_BatchSlot> and
                sizeof(CBCD_BatchSlot) == 5 * sizeof(std::uint32_t) + 2 * sizeof(std::uint16_t) and
                alignof(CBCD_BatchSlot) == alignof(std::uint32_t) and
                offsetof(CBCD_BatchSlot, block_end) == 0 and
                offsetof(CBCD_BatchSlot, first_entry) == sizeof(std::uint32_t) and
                offsetof(CBCD_BatchSlot, num_cells) == 2 * sizeof(std::uint32_t) and
                offsetof(CBCD_BatchSlot, num_incoming_faces) == 3 * sizeof(std::uint32_t) and
                offsetof(CBCD_BatchSlot, angle_set_id) == 4 * sizeof(std::uint32_t) and
                offsetof(CBCD_BatchSlot, block_size_x) == 5 * sizeof(std::uint32_t) and
                offsetof(CBCD_BatchSlot, num_stride_blocks) ==
                  5 * sizeof(std::uint32_t) + sizeof(std::uint16_t),
              "Device CBC batch-slot field layout changed.");
// NOLINTNEXTLINE(misc-redundant-expression): Profiles may share class extents deliberately.
static_assert(IsCBCDLaunchLayoutValid<CBCD_LAUNCH_CLASS_BYTES[0]>() and
                IsCBCDLaunchLayoutValid<CBCD_LAUNCH_CLASS_BYTES[1]>() and
                IsCBCDLaunchLayoutValid<CBCD_LAUNCH_CLASS_BYTES[2]>(),
              "Device CBC launch header, slots, entries or metadata extent changed.");

} // namespace opensn::gpu_kernel
