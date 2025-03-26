// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

/// Object holding a grouping.
class LBSGroup
{
public:
  int id;

public:
  LBSGroup() : id(-1) {}
  explicit LBSGroup(int id) : id(id) {}
};

} // namespace opensn
