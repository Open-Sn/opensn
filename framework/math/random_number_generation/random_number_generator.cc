#include "framework/math/random_number_generation/random_number_generator.h"

namespace opensn
{

RandomNumberGenerator::RandomNumberGenerator() : distribution_(0.0, 1.0)
{
  mt1993764_generator_.seed(0);
}

RandomNumberGenerator::RandomNumberGenerator(int seed) : distribution_(0.0, 1.0)
{
  mt1993764_generator_.seed(seed);
}

double
RandomNumberGenerator::Rand()
{
  return distribution_(mt1993764_generator_);
}

} // namespace opensn
