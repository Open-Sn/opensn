#include "framework/utils/chi_utils.h"
#include "framework/logging/chi_log_exceptions.h"
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>

#define scdouble static_cast<double>

namespace chi
{

std::string
StringLTrim(const std::string& s)
{
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

std::string
StringRTrim(const std::string& s)
{
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string
StringTrim(const std::string& s)
{
  return StringRTrim(StringLTrim(s));
}

std::vector<std::string>
StringSplit(const std::string& input, const std::string& delim)
{
  constexpr size_t NPOS = std::string::npos;
  std::vector<std::string> output;

  std::string remainder = input;
  size_t first_scope = remainder.find_first_of(delim);

  while (first_scope != NPOS)
  {
    if (first_scope != 0) output.push_back(remainder.substr(0, first_scope));

    remainder = remainder.substr(first_scope + delim.size(), NPOS);
    first_scope = remainder.find_first_of(delim);
  }
  output.push_back(remainder);

  return output;
}

std::string
StringUpToFirstReverse(const std::string& input, const std::string& search_string)
{
  constexpr size_t NPOS = std::string::npos;
  std::string output = input;
  const size_t last_scope = input.find_last_of(search_string);
  if (last_scope != NPOS) output = input.substr(last_scope + search_string.size(), NPOS);

  return output;
}

void
AssertReadibleFile(const std::string& file_name)
{
  std::ifstream file(file_name.c_str(), std::ifstream::in);
  ChiLogicalErrorIf(file.fail(),
                    "Failed to open file \"" + file_name +
                      "\"."
                      "Either the file does not exist or you do not have read permissions.");

  file.close();
}

std::string
PrintIterationProgress(const size_t current_iteration,
                       const size_t total_num_iterations,
                       const unsigned int num_intvls)
{
  typedef unsigned int uint;

  // Creating shorthand symbols for arguments
  const auto& i = current_iteration;
  const auto& I = num_intvls;
  const auto& N = total_num_iterations;

  // If on the first iteration then do nothing
  if (i == 0) return {};

  // Prepare an output stream
  std::stringstream output;
  output << std::fixed << std::setprecision(2) << std::setw(7);

  // If at the end, just print 100 and quit
  if ((i + 1) == N)
  {
    output << 100.0;
    return output.str();
  }

  const double dI = std::ceil(double(N) / I); // Interval size

  // std::modf is used to get the integral part
  // of a real value
  double x1;
  std::modf(double(i - 1) / dI, &x1);
  double x2;
  std::modf(double(i) / dI, &x2);

  if (uint(x2) != uint(x1))
  {
    output << x2 * (100.0 / I);
    return output.str();
  }

  return {};
}

std::vector<SubSetInfo>
MakeSubSets(size_t num_items, size_t desired_num_subsets)
{
  std::vector<SubSetInfo> ss_infos;
  const std::size_t div = std::floor(scdouble(num_items) / scdouble(desired_num_subsets));
  const std::size_t rem = num_items % desired_num_subsets;

  for (size_t i = 0; i < desired_num_subsets; ++i)
    ss_infos.push_back({0, 0, div});
  for (size_t j = 0; j < rem; ++j)
    ss_infos[j].ss_size += 1;

  size_t check_sum = 0;
  for (size_t i = 0; i < desired_num_subsets; ++i)
  {
    ss_infos[i].ss_begin = check_sum;
    check_sum += ss_infos[i].ss_size;
    ss_infos[i].ss_end = check_sum - 1;
  }

  return ss_infos;
}

} // namespace chi
