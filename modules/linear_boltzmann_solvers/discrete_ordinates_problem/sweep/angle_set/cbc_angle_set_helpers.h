// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace caribou
{
class Stream;
} // namespace caribou

namespace opensn
{

class CBC_AngleSet;

/// Get the caribou stream associated with the CBC_AngleSet
caribou::Stream& GetCBCAngleSetStream(CBC_AngleSet& angle_set);

/// Get the caribou stream associated with the CBC_AngleSet (const version)
const caribou::Stream& GetCBCAngleSetStream(const CBC_AngleSet& angle_set);

} // namespace opensn