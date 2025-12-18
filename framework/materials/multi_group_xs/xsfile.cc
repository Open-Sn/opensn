// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/multi_group_xs/xsfile.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

#include <cmath>
#include <unordered_set>

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
  bool read_num_groups = false;
  bool read_num_moments = false;
  bool read_num_precursors = false;
  std::unordered_set<std::string> seen_sections;

  const auto mark_unique_header = [&](const std::string& name, bool& flag)
  {
    OpenSnLogicalErrorIf(flag, "Duplicate top-level header \"" + name + "\" encountered.");
    flag = true;
  };
  const auto mark_unique_section = [&](const std::string& name)
  {
    OpenSnLogicalErrorIf(seen_sections.count(name) > 0,
                         "Duplicate section \"" + name + "_BEGIN\" encountered.");
    seen_sections.insert(name);
  };
  const auto require_num_groups = [&]()
  {
    OpenSnLogicalErrorIf(not read_num_groups,
                         "NUM_GROUPS must be specified before reading group-dependent data.");
  };
  const auto require_num_precursors = [&]()
  {
    OpenSnLogicalErrorIf(not read_num_precursors,
                         "NUM_PRECURSORS must be specified before reading precursor data.");
  };
  const auto IsFinite = [](double x) { return std::isfinite(x); };
  const auto IsFiniteVector = [&](const std::vector<double>& v)
  { return std::all_of(v.begin(), v.end(), IsFinite); };

  while (std::getline(file_, line))
  {
    ++line_number;
    std::istringstream line_stream(line);
    if (not(line_stream >> word))
      continue;
    if (word.front() == '#')
      continue;

    // Parse number of groups
    if (word == "NUM_GROUPS")
    {
      mark_unique_header("NUM_GROUPS", read_num_groups);
      int n_groups = 0;
      OpenSnLogicalErrorIf(not(static_cast<bool>(line_stream >> n_groups)),
                           "Failed parsing NUM_GROUPS value.");
      OpenSnLogicalErrorIf(n_groups <= 0, "The number of energy groups must be positive.");
      num_groups_ = n_groups;
      continue;
    }

    // Parse the number of scattering moments
    if (word == "NUM_MOMENTS")
    {
      mark_unique_header("NUM_MOMENTS", read_num_moments);
      int n_moments = 0;
      OpenSnLogicalErrorIf(not(static_cast<bool>(line_stream >> n_moments)),
                           "Failed parsing NUM_MOMENTS value.");
      OpenSnLogicalErrorIf(n_moments < 0, "The number of scattering moments must be non-negative.");
      scattering_order_ = std::max(0, n_moments - 1);
      continue;
    }

    // Parse the number of precursors species
    if (word == "NUM_PRECURSORS")
    {
      mark_unique_header("NUM_PRECURSORS", read_num_precursors);
      int n_prec = 0;
      OpenSnLogicalErrorIf(not(static_cast<bool>(line_stream >> n_prec)),
                           "Failed parsing NUM_PRECURSORS value.");
      OpenSnLogicalErrorIf(n_prec < 0,
                           "The number of delayed neutron precursors must be non-negative.");
      num_precursors_ = n_prec;
      precursors_.resize(num_precursors_);
      continue;
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
      {
        require_num_groups();
        mark_unique_section("GROUP_STRUCTURE");
        ReadGroupStructure("GROUP_STRUCTURE", e_bounds_, num_groups_, f, ls, ln);
      }
      if (fw == "INV_VELOCITY_BEGIN")
      {
        require_num_groups();
        mark_unique_section("INV_VELOCITY");
        Read1DData("INV_VELOCITY", inv_velocity_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsPositive(inv_velocity_),
                             "Only positive inverse velocity values are permitted.");
        OpenSnLogicalErrorIf(not IsFiniteVector(inv_velocity_),
                             "Inverse velocity values must be finite.");
      }
      if (fw == "VELOCITY_BEGIN" and inv_velocity_.empty())
      {
        require_num_groups();
        mark_unique_section("VELOCITY");
        Read1DData("VELOCITY", inv_velocity_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsPositive(inv_velocity_),
                             "Only positive velocity values are permitted.");
        OpenSnLogicalErrorIf(not IsFiniteVector(inv_velocity_), "Velocity values must be finite.");

        // Compute inverse velocity
        for (unsigned int g = 0; g < num_groups_; ++g)
          inv_velocity_[g] = 1.0 / inv_velocity_[g];
      }

      //
      // Read cross-section data
      //

      if (fw == "SIGMA_T_BEGIN")
      {
        require_num_groups();
        mark_unique_section("SIGMA_T");
        Read1DData("SIGMA_T", sigma_t_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(sigma_t_),
                             "Only non-negative total cross-section values are permitted.");
        OpenSnLogicalErrorIf(not IsFiniteVector(sigma_t_),
                             "Total cross section values must be finite.");
      } // if sigma_t

      if (fw == "SIGMA_A_BEGIN")
      {
        require_num_groups();
        mark_unique_section("SIGMA_A");
        Read1DData("SIGMA_A", sigma_a_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(sigma_a_),
                             "Only non-negative absorption cross-section values are permitted.");
        OpenSnLogicalErrorIf(not IsFiniteVector(sigma_a_),
                             "Absorption cross section values must be finite.");
      } // if sigma_a

      if (fw == "SIGMA_F_BEGIN")
      {
        require_num_groups();
        mark_unique_section("SIGMA_F");
        Read1DData("SIGMA_F", sigma_f_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(sigma_f_),
                             "Only non-negative fission cross-section values are permitted.");
        OpenSnLogicalErrorIf(not IsFiniteVector(sigma_f_),
                             "Fission cross section values must be finite.");
        if (not HasNonZero(sigma_f_))
        {
          log.Log0Warning() << "The fission cross section specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          sigma_f_.clear();
        }
      } // if sigma_f

      if (fw == "NU_SIGMA_F_BEGIN")
      {
        require_num_groups();
        mark_unique_section("NU_SIGMA_F");
        Read1DData("NU_SIGMA_F", nu_sigma_f_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(
          not IsNonNegative(nu_sigma_f_),
          "Only non-negative total fission multiplication cross-section values are permitted.");
        OpenSnLogicalErrorIf(not IsFiniteVector(nu_sigma_f_),
                             "Fission multiplication cross section values must be finite.");
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
        require_num_groups();
        mark_unique_section("NU");
        Read1DData("NU", nu_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(
          not std::all_of(nu_.begin(), nu_.end(), [](double x) { return x == 0.0 or x > 1.0; }),
          "Total fission neutron yield values must be either zero, or greater than one.");
        OpenSnLogicalErrorIf(not IsFiniteVector(nu_),
                             "Total fission neutron yield must be finite.");
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
        require_num_groups();
        mark_unique_section("NU_PROMPT");
        Read1DData("NU_PROMPT", nu_prompt_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not std::all_of(nu_prompt_.begin(),
                                             nu_prompt_.end(),
                                             [](double x) { return x == 0.0 or x > 1.0; }),
                             "Average prompt fission neutron yield values must be either zero, "
                             "or greater than one.");
        OpenSnLogicalErrorIf(not IsFiniteVector(nu_prompt_),
                             "Prompt fission neutron yield must be finite.");
        if (not HasNonZero(nu_prompt_))
        {
          log.Log0Warning() << "The prompt fission neutron yield specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          nu_prompt_.clear();
        }
      } // if nu_prompt

      if (fw == "NU_DELAYED_BEGIN")
      {
        require_num_groups();
        mark_unique_section("NU_DELAYED");
        Read1DData("NU_DELAYED", nu_delayed_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not IsNonNegative(nu_delayed_),
                             "Average delayed fission neutron yield values "
                             "must be non-negative.");
        OpenSnLogicalErrorIf(not IsFiniteVector(nu_delayed_),
                             "Delayed fission neutron yield must be finite.");
        if (not HasNonZero(nu_delayed_))
        {
          log.Log0Warning() << "The delayed fission neutron yield specified in "
                            << "\"" << file_name_ << "\" is uniformly zero... Clearing it.";
          nu_prompt_.clear();
        }
      } // if nu_delayed

      if (fw == "BETA_BEGIN")
      {
        require_num_groups();
        mark_unique_section("BETA");
        Read1DData("BETA", beta_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not std::all_of(beta_.begin(),
                                             beta_.end(),
                                             [](double x) { return x >= 0.0 and x <= 1.0; }),
                             "Delayed neutron fraction values must be in the range [0.0, 1.0].");
        OpenSnLogicalErrorIf(not IsFiniteVector(beta_),
                             "Delayed neutron fractions must be finite.");
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
        require_num_groups();
        mark_unique_section("CHI");
        Read1DData("CHI", chi_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(
          not HasNonZero(chi_),
          "The steady-state fission spectrum must have at least one non-zero value.");
        OpenSnLogicalErrorIf(not IsNonNegative(chi_),
                             "The steady-state fission spectrum must be non-negative.");
        OpenSnLogicalErrorIf(not IsFiniteVector(chi_),
                             "Steady-state fission spectrum values must be finite.");

        // Normalizing
        const auto sum = std::accumulate(chi_.begin(), chi_.end(), 0.0);
        std::transform(
          chi_.begin(), chi_.end(), chi_.begin(), [sum](double& x) { return x / sum; });
      } // if chi

      if (fw == "CHI_PROMPT_BEGIN")
      {
        require_num_groups();
        mark_unique_section("CHI_PROMPT");
        Read1DData("CHI_PROMPT", chi_prompt_, num_groups_, f, ls, ln);
        OpenSnLogicalErrorIf(not HasNonZero(chi_prompt_),
                             "The prompt fission spectrum must have at least one non-zero value.");
        OpenSnLogicalErrorIf(not IsNonNegative(chi_prompt_),
                             "The prompt fission spectrum must be non-negative.");
        OpenSnLogicalErrorIf(not IsFiniteVector(chi_prompt_),
                             "Prompt fission spectrum values must be finite.");

        // Normalizing
        const auto sum = std::accumulate(chi_prompt_.begin(), chi_prompt_.end(), 0.0);
        std::transform(chi_prompt_.begin(),
                       chi_prompt_.end(),
                       chi_prompt_.begin(),
                       [sum](double& x) { return x / sum; });

      } // if prompt chi

      if (num_precursors_ > 0 and fw == "CHI_DELAYED_BEGIN")
      {
        require_num_groups();
        require_num_precursors();
        mark_unique_section("CHI_DELAYED");
        // TODO: Should the be flipped to PRECURSOR_G_VAL?
        Read2DData("CHI_DELAYED",
                   "G_PRECURSOR_VAL",
                   emission_spectra_,
                   num_precursors_,
                   num_groups_,
                   f,
                   ls,
                   ln);

        for (unsigned int j = 0; j < num_precursors_; ++j)
        {
          OpenSnLogicalErrorIf(not HasNonZero(emission_spectra_[j]),
                               "Delayed emission spectrum for precursor " + std::to_string(j) +
                                 " must have at least one non-zero value.");
          OpenSnLogicalErrorIf(not IsNonNegative(emission_spectra_[j]),
                               "Delayed emission spectrum for precursor " + std::to_string(j) +
                                 " must be non-negative.");
          OpenSnLogicalErrorIf(not IsFiniteVector(emission_spectra_[j]),
                               "Delayed emission spectrum for precursor " + std::to_string(j) +
                                 " must be finite.");

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
          require_num_precursors();
          mark_unique_section("PRECURSOR_DECAY_CONSTANTS");
          Read1DData("PRECURSOR_DECAY_CONSTANTS", decay_constants_, num_precursors_, f, ls, ln);
          OpenSnLogicalErrorIf(not IsPositive(decay_constants_),
                               "Delayed neutron precursor decay constants must be positive.");
          OpenSnLogicalErrorIf(not IsFiniteVector(decay_constants_),
                               "Delayed neutron precursor decay constants must be finite.");
        } // if decay constants

        if (fw == "PRECURSOR_FRACTIONAL_YIELDS_BEGIN")
        {
          require_num_precursors();
          mark_unique_section("PRECURSOR_FRACTIONAL_YIELDS");
          Read1DData("PRECURSOR_FRACTIONAL_YIELDS", fractional_yields_, num_precursors_, f, ls, ln);
          OpenSnLogicalErrorIf(
            not HasNonZero(fractional_yields_),
            "Delayed neutron precursor fractional yields must contain at least one non-zero.");
          OpenSnLogicalErrorIf(
            not std::all_of(fractional_yields_.begin(),
                            fractional_yields_.end(),
                            [](double x) { return x >= 0.0 and x <= 1.0; }),
            "Delayed neutron precursor fractional yields must be in the range [0.0, 1.0].");
          OpenSnLogicalErrorIf(not IsFiniteVector(fractional_yields_),
                               "Delayed neutron precursor fractional yields must be finite.");

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
      {
        require_num_groups();
        mark_unique_section("TRANSFER_MOMENTS");
        ReadTransferMatrices(
          "TRANSFER_MOMENTS", transfer_matrices_, scattering_order_ + 1, num_groups_, f, ls, ln);
      }

      if (fw == "PRODUCTION_MATRIX_BEGIN")
      {
        require_num_groups();
        mark_unique_section("PRODUCTION_MATRIX");
        Read2DData("PRODUCTION_MATRIX",
                   "GPRIME_G_VAL",
                   production_matrix_,
                   num_groups_,
                   num_groups_,
                   f,
                   ls,
                   ln);
      }

      const bool recognized =
        fw == "GROUP_STRUCTURE_BEGIN" or fw == "INV_VELOCITY_BEGIN" or fw == "VELOCITY_BEGIN" or
        fw == "SIGMA_T_BEGIN" or fw == "SIGMA_A_BEGIN" or fw == "SIGMA_F_BEGIN" or
        fw == "NU_SIGMA_F_BEGIN" or fw == "NU_BEGIN" or fw == "NU_PROMPT_BEGIN" or
        fw == "NU_DELAYED_BEGIN" or fw == "BETA_BEGIN" or fw == "CHI_BEGIN" or
        fw == "CHI_PROMPT_BEGIN" or fw == "CHI_DELAYED_BEGIN" or
        fw == "PRECURSOR_DECAY_CONSTANTS_BEGIN" or fw == "PRECURSOR_FRACTIONAL_YIELDS_BEGIN" or
        fw == "TRANSFER_MOMENTS_BEGIN" or fw == "PRODUCTION_MATRIX_BEGIN";
      OpenSnLogicalErrorIf(not recognized and not fw.empty() and fw.front() != '#',
                           "Unknown top-level token \"" + fw + "\" encountered.");
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

  OpenSnLogicalErrorIf(not read_num_groups, "Cross-section file is missing required NUM_GROUPS.");
  OpenSnLogicalErrorIf(sigma_t_.empty() and sigma_a_.empty(),
                       "Cross-section file must provide at least one of SIGMA_T or SIGMA_A.");

  // Cross-consistency checks for fission data.
  if (not sigma_f_.empty() and not nu_.empty() and not nu_sigma_f_.empty())
  {
    constexpr double abs_tol = 1.0e-12;
    constexpr double rel_tol = 1.0e-8;
    for (size_t g = 0; g < num_groups_; ++g)
    {
      const double expected = nu_[g] * sigma_f_[g];
      const double diff = std::fabs(expected - nu_sigma_f_[g]);
      const double scale = std::max({1.0, std::fabs(expected), std::fabs(nu_sigma_f_[g])});
      OpenSnLogicalErrorIf(diff > abs_tol + rel_tol * scale,
                           "Inconsistent fission data at group " + std::to_string(g) +
                             ": NU_SIGMA_F differs from NU*SIGMA_F.");
    }
  }

  if (not nu_.empty() and not nu_prompt_.empty())
    OpenSnLogicalError("Ambiguous fission neutron yield. Both NU and NU_PROMPT were specified.");
  if (not chi_.empty() and not chi_prompt_.empty())
    OpenSnLogicalError("Ambiguous fission spectrum. Both CHI and CHI_PROMPT were specified.");
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
  OpenSnLogicalErrorIf(not std::getline(file, line),
                       "Unexpected EOF encountered while parsing " + keyword + ".");
  line_stream = std::istringstream(line);
  ++line_number;
  while (true)
  {
    std::string token;
    if (not(line_stream >> token))
    {
      OpenSnLogicalErrorIf(not std::getline(file, line),
                           "Unexpected EOF encountered while parsing " + keyword + ".");
      line_stream = std::istringstream(line);
      ++line_number;
      continue;
    }
    if (token == keyword + "_END")
      break;
    if (token.front() == '#')
    {
      OpenSnLogicalErrorIf(not std::getline(file, line),
                           "Unexpected EOF encountered while parsing " + keyword + ".");
      line_stream = std::istringstream(line);
      ++line_number;
      continue;
    }

    // Get data from current line
    std::istringstream token_stream(token);
    if (not(static_cast<bool>(token_stream >> bound)))
    {
      std::string msg("Unknown entry token \"");
      msg.append(token).append("\" encountered while parsing ").append(keyword).append(".");
      OpenSnLogicalError(msg);
    }
    OpenSnLogicalErrorIf(not std::isfinite(bound), "Group-structure values must be finite.");
    destination.push_back(bound);
    OpenSnLogicalErrorIf(count++ >= n_grps + 1,
                         "Too many entries encountered when parsing group structure.\n"
                         "The expected number of entries is " +
                           std::to_string(n_grps + 1) + ".");

    // Go to next line
    OpenSnLogicalErrorIf(not std::getline(file, line),
                         "Unexpected EOF encountered while parsing " + keyword + ".");
    line_stream = std::istringstream(line);
    ++line_number;
  }

  OpenSnLogicalErrorIf(destination.size() != n_grps + 1,
                       "Incorrect number of entries encountered when parsing group structure.\n"
                       "The expected number of entries is " +
                         std::to_string(n_grps + 1) + ", but " +
                         std::to_string(destination.size()) + " were provided.");
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
  std::vector<bool> seen(n_entries, false);

  std::string line;
  size_t count = 0;

  // Read the block
  OpenSnLogicalErrorIf(not std::getline(file, line),
                       "Unexpected EOF encountered while parsing " + keyword + ".");
  line_stream = std::istringstream(line);
  ++line_number;
  while (true)
  {
    std::string first_token;
    if (line_stream >> first_token)
    {
      if (first_token == keyword + "_END")
        break;

      if (first_token.front() != '#')
      {
        // Parse strict "<idx> <value>" format.
        std::size_t i = 0;
        double value = 0.0;
        try
        {
          i = static_cast<std::size_t>(std::stoul(first_token));
          OpenSnLogicalErrorIf(not(static_cast<bool>(line_stream >> value)),
                               "Failed parsing indexed 1D data line.");
        }
        catch (...)
        {
          std::string msg("Unknown entry token \"");
          msg.append(first_token)
            .append("\" encountered while parsing ")
            .append(keyword)
            .append(".");
          OpenSnLogicalError(msg);
        }

        OpenSnLogicalErrorIf(i >= n_entries,
                             "1D data index out of range. Found index " + std::to_string(i) +
                               " for expected size " + std::to_string(n_entries) + ".");
        OpenSnLogicalErrorIf(seen.at(i),
                             "Duplicate 1D data index encountered: " + std::to_string(i) + ".");
        OpenSnLogicalErrorIf(count >= n_entries,
                             "Too many entries encountered when parsing 1D data.\n"
                             "The expected number of entries is " +
                               std::to_string(n_entries) + ".");

        destination.at(i) = value;
        seen.at(i) = true;
        ++count;
      }
    }

    // Go to next line
    OpenSnLogicalErrorIf(not std::getline(file, line),
                         "Unexpected EOF encountered while parsing " + keyword + ".");
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
  std::unordered_set<std::size_t> seen;

  // Read the block
  OpenSnLogicalErrorIf(not std::getline(file, line),
                       "Unexpected EOF encountered while parsing " + keyword + ".");
  line_stream = std::istringstream(line);
  ++line_number;
  while (true)
  {
    if (not(line_stream >> word))
    {
      OpenSnLogicalErrorIf(not std::getline(file, line),
                           "Unexpected EOF encountered while parsing " + keyword + ".");
      line_stream = std::istringstream(line);
      ++line_number;
      continue;
    }

    if (word == keyword + "_END")
      break;

    // Check that this line contains an entry
    if (word == entry_prefix)
    {
      // Get data from current line
      OpenSnLogicalErrorIf(not(static_cast<bool>(line_stream >> i >> j >> value)),
                           "Failed parsing 2D data line for entry prefix \"" + entry_prefix +
                             "\".");
      OpenSnLogicalErrorIf(not std::isfinite(value), "2D data values must be finite.");
      if (entry_prefix == "G_PRECURSOR_VAL")
      {
        // File format is (group, precursor, value), destination storage is [precursor][group].
        OpenSnLogicalErrorIf(i >= n_cols,
                             "2D group index out of range. Found index " + std::to_string(i) +
                               " for expected size " + std::to_string(n_cols) + ".");
        OpenSnLogicalErrorIf(j >= n_rows,
                             "2D precursor index out of range. Found index " + std::to_string(j) +
                               " for expected size " + std::to_string(n_rows) + ".");
        const auto key = j * n_cols + i;
        OpenSnLogicalErrorIf(seen.count(key) > 0,
                             "Duplicate 2D data index encountered at (group=" + std::to_string(i) +
                               ", precursor=" + std::to_string(j) + ").");
        seen.insert(key);
        destination.at(j).at(i) = value; // hack
      }
      else
      {
        OpenSnLogicalErrorIf(i >= n_rows,
                             "2D row index out of range. Found index " + std::to_string(i) +
                               " for expected size " + std::to_string(n_rows) + ".");
        OpenSnLogicalErrorIf(j >= n_cols,
                             "2D column index out of range. Found index " + std::to_string(j) +
                               " for expected size " + std::to_string(n_cols) + ".");
        const auto key = i * n_cols + j;
        OpenSnLogicalErrorIf(seen.count(key) > 0,
                             "Duplicate 2D data index encountered at (" + std::to_string(i) + "," +
                               std::to_string(j) + ").");
        seen.insert(key);
        destination.at(i).at(j) = value;
      }
    }
    else if (word.front() != '#')
    {
      std::string msg("Unknown entry token \"");
      msg.append(word).append("\" encountered while parsing ").append(keyword).append(".");
      OpenSnLogicalError(msg);
    }

    // Go to next line
    OpenSnLogicalErrorIf(not std::getline(file, line),
                         "Unexpected EOF encountered while parsing " + keyword + ".");
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
  std::unordered_set<std::size_t> seen;

  // Read the block
  OpenSnLogicalErrorIf(not std::getline(file, line),
                       "Unexpected EOF encountered while parsing " + keyword + ".");
  line_stream = std::istringstream(line);
  ++line_number;
  while (true)
  {
    if (not(line_stream >> word))
    {
      OpenSnLogicalErrorIf(not std::getline(file, line),
                           "Unexpected EOF encountered while parsing " + keyword + ".");
      line_stream = std::istringstream(line);
      ++line_number;
      continue;
    }

    if (word == keyword + "_END")
      break;

    // Check that this line contains an entry
    if (word == "M_GFROM_GTO_VAL")
    {
      // Get data from current line
      OpenSnLogicalErrorIf(not(static_cast<bool>(line_stream >> ell >> gfrom >> gto >> value)),
                           "Failed parsing transfer-matrix data line.");
      OpenSnLogicalErrorIf(not std::isfinite(value), "Transfer-matrix values must be finite.");
      OpenSnLogicalErrorIf(ell >= n_moms,
                           "Transfer-matrix moment index out of range. Found index " +
                             std::to_string(ell) + " for expected size " + std::to_string(n_moms) +
                             ".");
      OpenSnLogicalErrorIf(gfrom >= n_grps,
                           "Transfer-matrix from-group index out of range. Found index " +
                             std::to_string(gfrom) + " for expected size " +
                             std::to_string(n_grps) + ".");
      OpenSnLogicalErrorIf(gto >= n_grps,
                           "Transfer-matrix to-group index out of range. Found index " +
                             std::to_string(gto) + " for expected size " + std::to_string(n_grps) +
                             ".");
      const auto key = (ell * n_grps + gfrom) * n_grps + gto;
      OpenSnLogicalErrorIf(
        seen.count(key) > 0,
        "Duplicate transfer-matrix entry encountered at (ell=" + std::to_string(ell) +
          ", gfrom=" + std::to_string(gfrom) + ", gto=" + std::to_string(gto) + ").");
      seen.insert(key);
      destination.at(ell).Insert(gto, gfrom, value);
    }
    else if (word.front() != '#')
    {
      std::string msg("Unknown entry token \"");
      msg.append(word).append("\" encountered while parsing ").append(keyword).append(".");
      OpenSnLogicalError(msg);
    }

    // Go to next line
    OpenSnLogicalErrorIf(not std::getline(file, line),
                         "Unexpected EOF encountered while parsing " + keyword + ".");
    line_stream = std::istringstream(line);
    ++line_number;
  }
}

} // namespace opensn
