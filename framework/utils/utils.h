#pragma once

#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

/**Miscellaneous utilities. These utilities should have no dependencies.*/
namespace opensn
{
/**Print the percentage completed based on the given interval.
 *
 * The function divides 100% into `num_intvls` intervals. If an
 * iteration passes an interval boundary then that interval percentage
 * will be printed.
 *
 * Specifying 10 intervals will print after each 10% is completed.
 * Specifying 4 intervals will print after each 25% is completed.*/
std::string PrintIterationProgress(size_t current_iteration,
                                   size_t total_num_iterations,
                                   unsigned int num_intvls = 10);

// #include <iostream>
// #include <string>

const std::string WHITESPACE = " \n\r\t\f\v";

/**Trims whitespace from the front of a string.*/
std::string StringLTrim(const std::string& s);

/**Trims whitespace from the back of a string.*/
std::string StringRTrim(const std::string& s);

/**Trims whitespace from the front and back of a string.*/
std::string StringTrim(const std::string& s);

/**Splits a string using the given delimiter. Consecutive delimiters
 * are treated as one.*/
std::vector<std::string> StringSplit(const std::string& input, const std::string& delim = " ");

/**The string portion, from the rear of the input string, up to
 * encountering the search_string.*/
std::string StringUpToFirstReverse(const std::string& input, const std::string& search_string);

void AssertReadibleFile(const std::string& file_name);

template <typename T, typename B>
bool
VectorListHas(const std::vector<T>& list, const B& val)
{
  return std::find(list.begin(), list.end(), val) != list.end();
}

struct SubSetInfo
{
  size_t ss_begin;
  size_t ss_end;
  size_t ss_size;
};
/**Subdivides a number of items (X) into a desired number of sub sets (Y). The
 * remainder of X/Y, i.e. r=X/Y obeys the indentity r < Y. These items will be
 * distributed to the first Y sub-sets. Example:
 * MakeSubSets(6659, 8) generates subsets of sizes
 * {833,833,833,832,832,832,832,832}.*/
std::vector<SubSetInfo> MakeSubSets(size_t num_items, size_t desired_num_subsets);

/**Popular and fast djb2a hashing algorithm.*/
inline constexpr uint32_t
hash_djb2a(const std::string_view sv)
{
  uint32_t hash{5381};
  for (unsigned char c : sv)
    hash = ((hash << 5) + hash) ^ c;

  return hash;
}

inline constexpr uint32_t operator"" _hash(const char* str, size_t len)
{
  return hash_djb2a(std::string_view{str, len});
}

template <typename T>
void
WriteBinaryValue(std::ofstream& output_file, T value)
{
  output_file.write((char*)&value, sizeof(T));
}
template <typename T>
T
ReadBinaryValue(std::ifstream& input_file)
{
  T value;
  input_file.read((char*)&value, sizeof(T));

  return value;
}
} // namespace opensn
