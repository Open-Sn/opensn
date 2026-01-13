// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

/// Object holding a grouping.
class LBSGroup
{
public:
  explicit LBSGroup(unsigned int id) : id(id) {}

  unsigned int id;
};

} // namespace opensn
