// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/xs.h"
#include <memory>

using namespace opensn;

namespace opensnlua
{

std::shared_ptr<MultiGroupXS>
XSCreate()
{
  auto xs = std::make_shared<MultiGroupXS>();
  opensn::multigroup_xs_stack.push_back(xs);
  const size_t index = opensn::multigroup_xs_stack.size() - 1;
  return xs;
}

std::shared_ptr<opensn::MultiGroupXS>
XSCreateSimpleOneGroup(double sigma_t, double c)
{
  auto xs = std::make_shared<MultiGroupXS>();
  xs->Initialize(sigma_t, c);
  return xs;
}

std::shared_ptr<opensn::MultiGroupXS>
XSLoadFromOpenSn(const std::string& file_name)
{
  auto xs = std::make_shared<MultiGroupXS>();
  xs->Initialize(file_name);
  return xs;
}

std::shared_ptr<opensn::MultiGroupXS>
XSLoadFromOpenMC(const std::string& file_name, const std::string& dataset_name, double temperature)
{
  auto xs = std::make_shared<MultiGroupXS>();
  xs->Initialize(file_name, dataset_name, temperature);
  return xs;
}

} // namespace opensnlua
