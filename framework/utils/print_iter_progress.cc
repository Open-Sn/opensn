#include "chi_utils.h"

#include <cmath>
#include <iomanip>
#include <sstream>

//###################################################################
/**Print the percentage completed based on the given interval.
 *
 * The function divides 100% into `num_intvls` intervals. If an
 * iteration passes an interval boundary then that interval percentage
 * will be printed.
 *
 * Specifying 10 intervals will print after each 10% is completed.
 * Specifying 4 intervals will print after each 25% is completed.*/
std::string chi::
  PrintIterationProgress(const size_t current_iteration,
                         const size_t total_num_iterations,
                         const unsigned int num_intvls/*=10*/)
{
  typedef unsigned int uint;

  // Creating shorthand symbols for arguments
  const auto& i = current_iteration;
  const auto& I = num_intvls;
  const auto& N = total_num_iterations;

  // If on the first iteration then do nothing
  if (i==0) return {};

  // Prepare an output stream
  std::stringstream output;
  output << std::fixed << std::setprecision(2) << std::setw(7);

  // If at the end, just print 100 and quit
  if ((i+1)==N) { output << 100.0; return output.str(); }

  const double dI = std::ceil(double(N)/I); //Interval size

  // std::modf is used to get the integral part
  // of a real value
  double x1; std::modf(double(i-1)/dI,&x1);
  double x2; std::modf(double(i  )/dI,&x2);

  if (uint(x2) != uint(x1)) { output << x2*(100.0/I); return output.str(); }

  return {};
}
