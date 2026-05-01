// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caliper/cali.h"

namespace opensn
{

class CaliperPhaseScope
{
public:
  CaliperPhaseScope(const char* name, int& depth) : name_(name), depth_(depth)
  {
    active_ = depth_++ == 0;
    if (active_)
      CALI_MARK_PHASE_BEGIN(name_);
  }

  ~CaliperPhaseScope()
  {
    --depth_;
    if (active_)
      CALI_MARK_PHASE_END(name_);
  }

  CaliperPhaseScope(const CaliperPhaseScope&) = delete;
  CaliperPhaseScope& operator=(const CaliperPhaseScope&) = delete;

private:
  const char* name_;
  int& depth_;
  bool active_ = false;
};

class CaliperRegionScope
{
public:
  CaliperRegionScope(const char* name, int& depth) : name_(name), depth_(depth)
  {
    active_ = depth_++ == 0;
    if (active_)
      CALI_MARK_BEGIN(name_);
  }

  ~CaliperRegionScope()
  {
    --depth_;
    if (active_)
      CALI_MARK_END(name_);
  }

  CaliperRegionScope(const CaliperRegionScope&) = delete;
  CaliperRegionScope& operator=(const CaliperRegionScope&) = delete;

private:
  const char* name_;
  int& depth_;
  bool active_ = false;
};

inline int&
CaliperSetupPhaseDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperSolvePhaseDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperPostPhaseDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperAGSScopeDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperWGSScopeDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperSteadyStateScopeDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperTransientScopeDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperPIScopeDepth()
{
  static int depth = 0;
  return depth;
}

inline int&
CaliperNLKEScopeDepth()
{
  static int depth = 0;
  return depth;
}

} // namespace opensn
