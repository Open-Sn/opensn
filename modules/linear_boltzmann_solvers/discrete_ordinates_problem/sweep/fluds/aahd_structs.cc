// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include <iomanip>

std::ostream&
operator<<(std::ostream& os, const opensn::AAHD_FaceNode& n)
{
  os << "AAHD_FaceNode(";
  os << "c=" << std::setw(5) << n.GetCellIndex() << ", ";
  os << "f=" << std::setw(2) << n.GetFaceIndex() << ", ";
  os << "fn=" << std::setw(2) << n.GetFaceNodeIndex();
  os << ")";
  return os;
}

std::ostream&
operator<<(std::ostream& out, const opensn::AAHD_NonLocalFaceNode& n)
{
  out << "NLEdgeNode(";
  const std::uint64_t& packed_face = std::get<2>(n.GetGlobalOrdering());
  out << std::get<0>(n.GetGlobalOrdering()) << "["
      << static_cast<uint16_t>((packed_face >> (3 * 16)) & 0xFFFF) << ", "
      << static_cast<uint16_t>((packed_face >> (2 * 16)) & 0xFFFF) << "]";
  out << " <-> ";
  out << std::get<1>(n.GetGlobalOrdering()) << "["
      << static_cast<uint16_t>((packed_face >> (1 * 16)) & 0xFFFF) << ", "
      << static_cast<uint16_t>((packed_face >> (0 * 16)) & 0xFFFF) << "]";
  out << ")";
  return out;
}

std::ostream&
operator<<(std::ostream& os, const opensn::AAHD_NodeIndex& e)
{
  os << "AAHD_NodeIndex(";
  if (e.IsUndefined())
  {
    os << "undefined";
  }
  else
  {
    os << "index=" << std::setw(6) << e.GetIndex() << ", ";
    if (e.IsOutGoing())
      os << "outgoing, ";
    else
      os << "incoming, ";
    if (e.IsBoundary())
    {
      os << "is_boundary=true";
    }
    else
    {
      if (e.IsDelayed())
        os << "is_delayed= true, ";
      else
        os << "is_delayed=false, ";
      if (e.IsLocal())
        os << "is_local=true";
      else
        os << "loc=" << e.GetRank();
    }
  }
  os << ")";
  return os;
}
