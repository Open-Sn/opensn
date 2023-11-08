#pragma once

#include <random>

namespace chi_math
{
// #########################################################
/**Random number generator based on threefry.*/
class RandomNumberGenerator
{
private:
  std::mt19937_64 mt1993764_generator_;
  std::uniform_real_distribution<double> distribution_;

public:
  /**Default constructor. Seeds the generator with a zero.*/
  RandomNumberGenerator();
  /**Constructor where a seed is supplied.*/
  RandomNumberGenerator(int seed);
  /**Generates a random number with the default distribution.*/
  double Rand();
};
} // namespace chi_math
