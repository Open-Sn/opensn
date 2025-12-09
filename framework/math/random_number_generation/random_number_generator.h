// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <random>

namespace opensn
{

/// Random number generator based on threefry
class RandomNumberGenerator
{
public:
  /// Default constructor. Seeds the generator with a zero.
  RandomNumberGenerator();

  /// Constructor where a seed is supplied.
  RandomNumberGenerator(int seed);

  /// Generates a random number with the default distribution.
  double Rand();

private:
  std::mt19937_64 mt1993764_generator_;
  std::uniform_real_distribution<double> distribution_;
};

} // namespace opensn
