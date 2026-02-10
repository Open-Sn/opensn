// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/xsfile.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

XSFile::XSFile(const std::string& file_name)
  : file_name_(file_name),
    file_(file_name),
    num_groups_(0),
    scattering_order_(0),
    num_precursors_(0)
{
  OpenSnLogicalErrorIf(not file_.is_open(),
                       "Failed to open cross-section file " + file_name_ + ".");
  opensn::log.Log() << "Reading OpenSn cross-section file \"" << file_name_ << "\"\n";
}

void
XSFile::Read()
{
  std::string word, line;
  size_t line_number = 0;
  while (std::getline(file_, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    // Parse number of groups
    if (word == "NUM_GROUPS")
    {
      int n_groups = 0;
      line_stream >> n_groups;
      OpenSnLogicalErrorIf(n_groups <= 0, "The number of energy groups must be positive.");
      num_groups_ = n_groups;
    }

    // Parse the number of scattering moments
    if (word == "NUM_MOMENTS")
    {
      int n_moments = 0;
      line_stream >> n_moments;
      OpenSnLogicalErrorIf(n_moments < 0, "The number of scattering moments must be non-negative.");
      scattering_order_ = std::max(0, n_moments - 1);
    }

    // Parse the number of precursors species
    if (word == "NUM_PRECURSORS")
    {
      int n_prec = 0;
      line_stream >> n_prec;
      OpenSnLogicalErrorIf(n_prec < 0,
                           "The number of delayed neutron precursors must be non-negative.");
      num_precursors_ = n_prec;
      precursors_.resize(num_precursors_);
    }

    // Parse nuclear data
    try
    {
      auto& ln = line_number;
      auto& ls = line_stream;
      auto& f = file_;
      auto& fw = word;

      //
      // Group Structure Data
      //

      if (fw == "GROUP_STRUCTURE_BEGIN")
        ReadGroupStructure("GROUP_STRUCTURE", e_bounds_, num_groups_, f, ls, ln);
      if (fw == "INV_VELOCITY_BEGIN")
      {
        Read1DData("INV_VELOCITY", inv_velocity_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsPositive(inv_velocity_),
                             "Only positive inverse velocity values are permitted.");
      }
      if (fw == "VELOCITY_BEGIN" and inv_velocity_.empty())
      {
        Read1DData("VELOCITY", inv_velocity_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsPositive(inv_velocity_),
                             "Only positive velocity values are permitted.");

        // Compute inverse velocity
        for (unsigned int g = 0; g < num_groups_; ++g)
          inv_velocity_[g] = 1.0 / inv_velocity_[g];
      }

      //
      // Read cross-section data
      //

      if (fw == "SIGMA_T_BEGIN")
      {
        Read1DData("SIGMA_T", sigma_t_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(sigma_t_),
                             "Only non-negative total cross-section values are permitted.");
      } // if sigma_t

      if (fw == "SIGMA_A_BEGIN")
      {
        Read1DData("SIGMA_A", sigma_a_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(sigma_a_),
                             "Only non-negative absorption cross-section values are permitted.");
      } // if sigma_a

      if (fw == "SIGMA_F_BEGIN")
      {
        Read1DData("SIGMA_F", sigma_f_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(sigma_f_),
                             "Only non-negative fission cross-section values are permitted.");
        if (not HasNonZero(sigma_f_))
        {
          log.Log0Warning() << "The fission cross section specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          sigma_f_.clear();
        }
      } // if sigma_f

      if (fw == "NU_SIGMA_F_BEGIN")
      {
        Read1DData("NU_SIGMA_F", nu_sigma_f_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(
          not IsNonNegative(nu_sigma_f_),
          "Only non-negative total fission multiplication cross-section values are permitted.");
        if (not HasNonZero(nu_sigma_f_))
        {
          log.Log0Warning() << "The production cross section specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          nu_sigma_f_.clear();
        }
      } // if nu_sigma_f

      //
      // Read neutrons per fission data
      //

      if (fw == "NU_BEGIN")
      {
        Read1DData("NU", nu_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(
          not std::all_of(nu_.begin(), nu_.end(), [](double x) { return x == 0.0 or x > 1.0; }),
          "Total fission neutron yield values must be either zero, or greater than one.");
        if (not HasNonZero(nu_))
        {
          log.Log0Warning() << "The total fission neutron yield specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          nu_.clear();
        }

        // Compute prompt/delayed nu, if needed
        if (num_precursors_ > 0 and not nu_.empty() and not beta_.empty() and nu_prompt_.empty() and
            nu_delayed_.empty())
        {
          nu_prompt_.assign(num_groups_, 0.0);
          nu_delayed_.assign(num_groups_, 0.0);
          for (unsigned int g = 0; g < num_groups_; ++g)
          {
            nu_prompt_[g] = (1.0 - beta_[g]) * nu_[g];
            nu_delayed_[g] = beta_[g] * nu_[g];
          }
        }
      } // if nu

      if (fw == "NU_PROMPT_BEGIN")
      {
        Read1DData("NU_PROMPT", nu_prompt_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not std::all_of(nu_prompt_.begin(),
                                             nu_prompt_.end(),
                                             [](double x) { return x == 0.0 or x > 1.0; }),
                             "Average prompt fission neutron yield values must be either zero, "
                             "or greater than one.");
        if (not HasNonZero(nu_prompt_))
        {
          log.Log0Warning() << "The prompt fission neutron yield specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          nu_prompt_.clear();
        }
      } // if nu_prompt

      if (fw == "NU_DELAYED_BEGIN")
      {
        Read1DData("NU_DELAYED", nu_delayed_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(nu_delayed_),
                             "Average delayed fission neutron yield values "
                             "must be non-negative.");
        if (not HasNonZero(nu_delayed_))
        {
          log.Log0Warning() << "The delayed fission neutron yield specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          nu_prompt_.clear();
        }
      } // if nu_delayed

      if (fw == "BETA_BEGIN")
      {
        Read1DData("BETA", beta_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not std::all_of(beta_.begin(),
                                             beta_.end(),
                                             [](double x) { return x >= 0.0 and x <= 1.0; }),
                             "Delayed neutron fraction values must be in the range [0.0, 1.0].");
        if (not HasNonZero(beta_))
        {
          log.Log0Warning() << "The delayed neutron fraction specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          beta_.clear();
        }

        // Compute prompt/delayed nu, if needed
        if (num_precursors_ > 0 and not nu_.empty() and not beta_.empty() and nu_prompt_.empty() and
            nu_delayed_.empty())
        {
          nu_prompt_.assign(num_groups_, 0.0);
          nu_delayed_.assign(num_groups_, 0.0);
          for (unsigned int g = 0; g < num_groups_; ++g)
          {
            nu_prompt_[g] = (1.0 - beta_[g]) * nu_[g];
            nu_delayed_[g] = beta_[g] * nu_[g];
          }
        }
      } // if beta

      //
      // Read fission/emission spectra
      //

      if (fw == "CHI_BEGIN")
      {
        Read1DData("CHI", chi_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(
          not HasNonZero(chi_),
          "The steady-state fission spectrum must have at least one non-zero value.");
        OpenSnLogicalErrorIf(not IsNonNegative(chi_),
                             "The steady-state fission spectrum must be non-negative.");

        // Normalizing
        const auto sum = std::accumulate(chi_.begin(), chi_.end(), 0.0);
        std::transform(
          chi_.begin(), chi_.end(), chi_.begin(), [sum](double& x) { return x / sum; });
      } // if chi

      if (fw == "CHI_PROMPT_BEGIN")
      {
        Read1DData("CHI_PROMPT", chi_prompt_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not HasNonZero(chi_prompt_),
                             "The prompt fission spectrum must have at least one non-zero value.");
        OpenSnLogicalErrorIf(not IsNonNegative(chi_prompt_),
                             "The prompt fission spectrum must be non-negative.");

        // Normalizing
        const auto sum = std::accumulate(chi_prompt_.begin(), chi_prompt_.end(), 0.0);
        std::transform(chi_prompt_.begin(),
                       chi_prompt_.end(),
                       chi_prompt_.begin(),
                       [sum](double& x) { return x / sum; });

      } // if prompt chi

      if (num_precursors_ > 0 and fw == "CHI_DELAYED_BEGIN")
      {
        // TODO: Should the be flipped to PRECURSOR_G_VAL?
        Read2DData("CHI_DELAYED",
                   "G_PRECURSOR_VAL",
                   emission_spectra_,
                   num_precursors_,
                   num_groups_,
                   f,
                   ls,
                   ln);

        for (size_t j = 0; j < num_precursors_; ++j)
        {
          OpenSnLogicalErrorIf(not HasNonZero(emission_spectra_[j]),
                               "Delayed emission spectrum for precursor " + std::to_string(j) +
                                 " must have at least one non-zero value.");
          OpenSnLogicalErrorIf(not IsNonNegative(emission_spectra_[j]),
                               "Delayed emission spectrum for precursor " + std::to_string(j) +
                                 " must be non-negative.");

          // normalizing
          const auto sum =
            std::accumulate(emission_spectra_[j].begin(), emission_spectra_[j].end(), 0.0);
          std::transform(emission_spectra_[j].begin(),
                         emission_spectra_[j].end(),
                         emission_spectra_[j].begin(),
                         [sum](double& x) { return x / sum; });
        }
      } // if delayed chi

      //
      // Read delayed neutron precursor data
      //

      if (num_precursors_ > 0)
      {
        if (fw == "PRECURSOR_DECAY_CONSTANTS_BEGIN")
        {
          Read1DData("PRECURSOR_DECAY_CONSTANTS", decay_constants_, num_precursors_, f, ls, ln);
          OpenSnLogicalErrorIf(not IsPositive(decay_constants_),
                               "Delayed neutron precursor decay constants must be positive.");
        } // if decay constants

        if (fw == "PRECURSOR_FRACTIONAL_YIELDS_BEGIN")
        {
          Read1DData("PRECURSOR_FRACTIONAL_YIELDS", fractional_yields_, num_precursors_, f, ls, ln);
          OpenSnLogicalErrorIf(
            not HasNonZero(fractional_yields_),
            "Delayed neutron precursor fractional yields must contain at least one non-zero.");
          OpenSnLogicalErrorIf(
            not std::all_of(fractional_yields_.begin(),
                            fractional_yields_.end(),
                            [](double x) { return x >= 0.0 and x <= 1.0; }),
            "Delayed neutron precursor fractional yields must be in the range [0.0, 1.0].");

          // Normalizing
          const auto sum =
            std::accumulate(fractional_yields_.begin(), fractional_yields_.end(), 0.0);
          std::transform(fractional_yields_.begin(),
                         fractional_yields_.end(),
                         fractional_yields_.begin(),
                         [sum](double& x) { return x / sum; });
        } // if precursor yield
      }

      //
      // Read transfer data
      //

      if (fw == "TRANSFER_MOMENTS_BEGIN")
        ReadTransferMatrices(
          "TRANSFER_MOMENTS", transfer_matrices_, scattering_order_ + 1, num_groups_, f, ls, ln);

      if (fw == "PRODUCTION_MATRIX_BEGIN")
        Read2DData("PRODUCTION_MATRIX",
                   "GPRIME_G_VAL",
                   production_matrix_,
                   num_groups_,
                   num_groups_,
                   f,
                   ls,
                   ln);
    } // try

    catch (const std::runtime_error& err)
    {
      throw std::runtime_error("Error reading OpenSn cross-section file "
                               "\"" +
                               file_name_ + "\".\n" + "Line number " + std::to_string(line_number) +
                               "\n" + err.what());
    }

    catch (const std::logic_error& err)
    {
      throw std::logic_error("Error reading OpenSn cross-section file "
                             "\"" +
                             file_name_ + "\".\n" + "Line number " + std::to_string(line_number) +
                             "\n" + err.what());
    }

    catch (...)
    {
      throw std::runtime_error("Unknown error encountered.");
    }

    word = "";
  } // while not EOF, read each lines
}

void
XSFile::ReadGroupStructure(const std::string& keyword,
                           std::vector<double>& destination,
                           const unsigned int n_grps,
                           std::ifstream& file,
                           std::istringstream& line_stream,
                           size_t& line_number)
{
  destination.reserve(n_grps + 1);

  std::string line;
  double bound = 0.0;
  size_t count = 0;

  // Read the block
  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;
  while (line != keyword + "_END")
  {
    // Get data from current line
    line_stream >> bound;
    destination.push_back(bound);
    OpenSnLogicalErrorIf(count++ >= n_grps + 1,
                         "Too many entries encountered when parsing group structure.\n"
                         "The expected number of entries is " +
                           std::to_string(n_grps + 1) + ".");

    // Go to next line
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}

void
XSFile::Read1DData(const std::string& keyword,
                   std::vector<double>& destination,
                   const std::size_t n_entries,
                   std::ifstream& file,
                   std::istringstream& line_stream,
                   size_t& line_number)
{
  destination.assign(n_entries, 0.0);

  std::string line;
  int i = 0;
  double value = 0.0;
  size_t count = 0;

  // Read the bloc
  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;
  while (line != keyword + "_END")
  {
    // get data from current line
    line_stream >> i >> value;
    destination.at(i) = value;
    OpenSnLogicalErrorIf(count++ >= n_entries,
                         "Too many entries encountered when parsing 1D data.\n"
                         "The expected number of entries is " +
                           std::to_string(n_entries) + ".");

    // Go to next line
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}

void
XSFile::Read2DData(const std::string& keyword,
                   const std::string& entry_prefix,
                   std::vector<std::vector<double>>& destination,
                   const size_t n_rows,
                   const size_t n_cols,
                   std::ifstream& file,
                   std::istringstream& line_stream,
                   size_t& line_number)
{
  destination.assign(n_rows, std::vector<double>(n_cols, 0.0));

  std::string word, line;
  double value = 0.0;
  size_t i = 0, j = 0;

  // Read the block
  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;
  while (line != keyword + "_END")
  {
    // Check that this line contains an entry
    line_stream >> word;
    if (word == entry_prefix)
    {
      // Get data from current line
      line_stream >> i >> j >> value;
      if (entry_prefix == "G_PRECURSOR_VAL")
        destination.at(j).at(i) = value; // hack
      else
        destination.at(i).at(j) = value;
    }

    // Go to next line
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}

// Lambda for reading transfer matrix data.
void
XSFile::ReadTransferMatrices(const std::string& keyword,
                             std::vector<SparseMatrix>& destination,
                             const size_t n_moms,
                             const size_t n_grps,
                             std::ifstream& file,
                             std::istringstream& line_stream,
                             size_t& line_number)
{
  destination.assign(n_moms, SparseMatrix(n_grps, n_grps));

  std::string word, line;
  double value = 0.0;
  size_t ell = 0, gto = 0, gfrom = 0;

  // Read the block
  std::getline(file, line);
  line_stream = std::istringstream(line);
  ++line_number;
  while (line != keyword + "_END")
  {
    // Check that this line contains an entry
    line_stream >> word;
    if (word == "M_GFROM_GTO_VAL")
    {
      // Get data from current line
      line_stream >> ell >> gfrom >> gto >> value;
      destination.at(ell).Insert(gto, gfrom, value);
    }

    // Go to next line
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
  }
}

} // namespace opensn
