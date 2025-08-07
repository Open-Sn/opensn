// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include <vector>
#include <string>

namespace opensn
{

class XSFile
{
public:
  /**
   * Create a cross-section file
   *
   * @param file_name File name
   */
  XSFile(const std::string& file_name);

  /// Read the file
  MultiGroupXS Read();

private:
  /// Reading group structure data.
  void ReadGroupStructure(const std::string& keyword,
                          std::vector<double>& destination,
                          const size_t n_grps,
                          std::ifstream& file,
                          std::istringstream& line_stream,
                          size_t& line_number);

  /// Lambda for reading vector data.
  void Read1DData(const std::string& keyword,
                  std::vector<double>& destination,
                  const size_t n_entries,
                  std::ifstream& file,
                  std::istringstream& line_stream,
                  size_t& line_number);

  /// Lambda for reading 2D data
  void Read2DData(const std::string& keyword,
                  const std::string& entry_prefix,
                  std::vector<std::vector<double>>& destination,
                  const size_t n_rows,
                  const size_t n_cols,
                  std::ifstream& file,
                  std::istringstream& line_stream,
                  size_t& line_number);

  /// Lambda for reading transfer matrix data.
  void ReadTransferMatrices(const std::string& keyword,
                            std::vector<SparseMatrix>& destination,
                            const size_t n_moms,
                            const size_t n_grps,
                            std::ifstream& file,
                            std::istringstream& line_stream,
                            size_t& line_number);

  std::string file_name_;
};

} // namespace opensn
