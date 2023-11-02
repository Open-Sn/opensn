#include "chi_utils.h"

#include <fstream>

#include "chi_log_exceptions.h"

namespace chi
{
// #################################################################
std::string
StringLTrim(const std::string& s)
{
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

// #################################################################
std::string
StringRTrim(const std::string& s)
{
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

// #################################################################
std::string
StringTrim(const std::string& s)
{
  return StringRTrim(StringLTrim(s));
}

// #################################################################
std::vector<std::string>
StringSplit(const std::string& input, const std::string& delim /*=" "*/)
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

// #################################################################
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
} // namespace chi