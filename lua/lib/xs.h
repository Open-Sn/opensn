// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/multi_group_xs/multi_group_xs.h"

namespace opensnlua
{

std::shared_ptr<opensn::MultiGroupXS> XSCreate();

std::shared_ptr<opensn::MultiGroupXS> XSCreateSimpleOneGroup(double sigma_t, double c);

std::shared_ptr<opensn::MultiGroupXS> XSLoadFromOpenSn(const std::string& file_name);

std::shared_ptr<opensn::MultiGroupXS>
XSLoadFromOpenMC(const std::string& file_name, const std::string& dataset_name, double temperature);

} // namespace opensnlua
