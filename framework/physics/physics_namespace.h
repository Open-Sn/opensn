// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

enum class OperationType
{
  SINGLE_VALUE = 0,
  FROM_ARRAY = 1,
  SIMPLE_ONE_GROUP = 20,
  EXISTING = 22,
  OPENSN_XSFILE = 23
};

} // namespace opensn
