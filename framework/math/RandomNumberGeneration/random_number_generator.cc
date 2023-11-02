#include "random_number_generator.h"

/**Default constructor. Seeds the generator with a zero.*/
chi_math::RandomNumberGenerator::RandomNumberGenerator() : distribution_(0.0, 1.0)
{
  mt1993764_generator_.seed(0);
}

/**Constructor where a seed is supplied.*/
chi_math::RandomNumberGenerator::RandomNumberGenerator(int seed) : distribution_(0.0, 1.0)
{
  mt1993764_generator_.seed(seed);
}

/**Generates a random number with the default distribution.*/
double
chi_math::RandomNumberGenerator::Rand()
{
  return distribution_(mt1993764_generator_);
}