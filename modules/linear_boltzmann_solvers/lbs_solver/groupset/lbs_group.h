// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{
namespace lbs
{
/**Object holding a grouping.*/
class LBSGroup
{
public:
  int id_;

public:
  LBSGroup() : id_(-1) {}
  explicit LBSGroup(int id) : id_(id) {}
};

} // namespace lbs
} // namespace opensn
