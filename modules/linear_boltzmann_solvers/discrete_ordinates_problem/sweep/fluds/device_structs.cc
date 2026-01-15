// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/device_structs.h"
#include <iomanip>

std::ostream&
std::operator<<(std::ostream& out, const opensn::DeviceFaceNode& n)
{
  out << "DeviceFaceNode(";
  out << "c=" << std::setw(5) << n.GetCellIndex() << ", ";
  out << "f=" << std::setw(2) << n.GetFaceIndex() << ", ";
  out << "fn=" << std::setw(2) << n.GetFaceNodeIndex();
  out << ")";
  return out;
}

std::ostream&
std::operator<<(std::ostream& out, const opensn::DeviceNodeIndex& e)
{
  out << "DeviceNodeIndex(";
  if (e.IsUndefined())
  {
    out << "undefined";
  }
  else
  {
    out << "index=" << std::setw(6) << e.GetIndex() << ", ";
    if (e.IsOutgoing())
      out << "outgoing, ";
    else
      out << "incoming, ";
    if (e.IsDelayed())
      out << "is_delayed=true, ";
    else
      out << "is_delayed=false, ";
    if (e.IsBoundary())
      out << "is_boundary=true, ";
    else
      out << "is_boundary=false, ";
    if (e.IsLocal())
      out << "is_local=true";
    else
      out << "is_local=false";
  }
  out << ")";
  return out;
} // namespace std