// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include <vector>
#include <string>
#include <fstream>

namespace opensn
{

/**
 * Class for reading cross-section data from a native OpenSn file.
 */
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
  void Read();

  std::string file_name_;
  std::ifstream file_;
  size_t num_groups_;
  unsigned int scattering_order_;
  size_t num_precursors_;
  std::vector<MultiGroupXS::Precursor> precursors_;
  std::vector<double> inv_velocity_;
  std::vector<double> e_bounds_;
  std::vector<double> sigma_t_;
  std::vector<double> sigma_a_;
  std::vector<double> sigma_f_;
  std::vector<double> nu_;
  std::vector<double> nu_prompt_;
  std::vector<double> nu_delayed_;
  std::vector<double> beta_;
  std::vector<double> chi_prompt_;
  std::vector<double> nu_sigma_f_;
  std::vector<double> chi_;
  std::vector<double> decay_constants_;
  std::vector<double> fractional_yields_;
  std::vector<std::vector<double>> emission_spectra_;
  std::vector<SparseMatrix> transfer_matrices_;
  std::vector<std::vector<double>> production_matrix_;

private:
  /// Reading group structure data.
  void ReadGroupStructure(const std::string& keyword,
                          std::vector<double>& destination,
                          size_t n_grps,
                          std::ifstream& file,
                          std::istringstream& line_stream,
                          size_t& line_number);

  /// Lambda for reading vector data.
  void Read1DData(const std::string& keyword,
                  std::vector<double>& destination,
                  size_t n_entries,
                  std::ifstream& file,
                  std::istringstream& line_stream,
                  size_t& line_number);

  /// Lambda for reading 2D data
  void Read2DData(const std::string& keyword,
                  const std::string& entry_prefix,
                  std::vector<std::vector<double>>& destination,
                  size_t n_rows,
                  size_t n_cols,
                  std::ifstream& file,
                  std::istringstream& line_stream,
                  size_t& line_number);

  /// Lambda for reading transfer matrix data.
  void ReadTransferMatrices(const std::string& keyword,
                            std::vector<SparseMatrix>& destination,
                            size_t n_moms,
                            size_t n_grps,
                            std::ifstream& file,
                            std::istringstream& line_stream,
                            size_t& line_number);
};

} // namespace opensn
