// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <vector>
#include <cstddef>

namespace opensn
{

class LBSGroupset;

enum SourceType
{
  APPLY_FIXED_SOURCES = (1 << 0),
  APPLY_WGS_SCATTER_SOURCES = (1 << 1),
  APPLY_AGS_SCATTER_SOURCES = (1 << 2),
  APPLY_WGS_FISSION_SOURCES = (1 << 3),
  APPLY_AGS_FISSION_SOURCES = (1 << 4),
  SUPPRESS_WG_SCATTER = (1 << 5),
  ZERO_INCOMING_DELAYED_PSI = (1 << 6)
};

/// SourceFlags is a combination of `SourceType`s
struct SourceFlags
{
  SourceFlags() : flags_(0) {}
  SourceFlags(SourceType type) : flags_(type) {}

  SourceFlags& operator|=(SourceType type)
  {
    flags_ |= type;
    return *this;
  }

  SourceFlags& operator|=(const SourceFlags& src)
  {
    flags_ |= src.flags_;
    return *this;
  }

  bool Empty() const { return flags_ == 0; }

  void Unset(SourceType type) { flags_ &= ~type; }

  bool operator&(const SourceType& type) const { return flags_ & type; }

private:
  int flags_;
};

inline SourceFlags
operator|(SourceFlags s1, SourceFlags s2)
{
  SourceFlags src = s1;
  src |= s2;
  return src;
}

inline SourceFlags
operator|(SourceType f1, SourceType f2)
{
  SourceFlags src;
  src |= f1;
  src |= f2;
  return src;
}

using SetSourceFunction = std::function<void(const LBSGroupset& groupset,
                                             std::vector<double>& q,
                                             const std::vector<double>& phi,
                                             const SourceFlags source_flags)>;

} // namespace opensn
