// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/random_number_generation/random_number_generator.h"

namespace opensn
{
namespace
{
std::random_device rd;
}

RandomNumberGenerator::RandomNumberGenerator() : mt1993764_generator_(rd()), distribution_(0.0, 1.0)
{
}

RandomNumberGenerator::RandomNumberGenerator(int seed)
  : mt1993764_generator_(seed), distribution_(0.0, 1.0)
{
}

double
RandomNumberGenerator::Rand()
{
  return distribution_(mt1993764_generator_);
}

} // namespace opensn
