// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include <iomanip>

std::ostream&
std::operator<<(std::ostream& out, const opensn::AAHD_FaceNode& n)
{
  out << "AAHD_FaceNode(";
  out << "c=" << std::setw(5) << n.GetCellIndex() << ", ";
  out << "f=" << std::setw(2) << n.GetFaceIndex() << ", ";
  out << "fn=" << std::setw(2) << n.GetFaceNodeIndex();
  out << ")";
  return out;
}

std::ostream&
std::operator<<(std::ostream& out, const opensn::AAHD_NonLocalFaceNode& n)
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
std::operator<<(std::ostream& out, const opensn::AAHD_NodeIndex& e)
{
  out << "AAHD_NodeIndex(";
  if (e.IsUndefined())
  {
    out << "undefined";
  }
  else
  {
    out << "index=" << std::setw(6) << e.GetIndex() << ", ";
    if (e.IsOutGoing())
      out << "outgoing, ";
    else
      out << "incoming, ";
    if (e.IsBoundary())
    {
      out << "is_boundary=true";
    }
    else
    {
      if (e.IsDelayed())
        out << "is_delayed= true, ";
      else
        out << "is_delayed=false, ";
      if (e.IsLocal())
        out << "is_local=true";
      else
        out << "is_local=false";
    }
  }
  out << ")";
  return out;
}
