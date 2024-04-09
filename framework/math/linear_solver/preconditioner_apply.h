// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

template <class PCType, class VecType>
int PreconditionerApplication(PCType pc, VecType vector, VecType action);

} // namespace opensn
