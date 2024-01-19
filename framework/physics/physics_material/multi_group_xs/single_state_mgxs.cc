#include "framework/physics/physics_material/multi_group_xs/single_state_mgxs.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <string>
#include <numeric>

namespace opensn
{

void
SingleStateMGXS::Clear()
{
  num_groups_ = 0;
  scattering_order_ = 0;
  num_precursors_ = 0;
  is_fissionable_ = false;

  sigma_t_.clear();
  sigma_a_.clear();
  transfer_matrices_.clear();

  sigma_f_.clear();
  nu_sigma_f_.clear();
  nu_prompt_sigma_f_.clear();
  nu_delayed_sigma_f_.clear();
  production_matrix_.clear();
  precursors_.clear();

  inv_velocity_.clear();

  // Diffusion quantities
  diffusion_initialized_ = false;
  sigma_tr_.clear();
  diffusion_coeff_.clear();
  sigma_r_.clear();
  sigma_s_gtog_.clear();

  // Monte-Carlo quantities
  scattering_initialized_ = false;
  cdf_gprime_g_.clear();
  scat_angles_gprime_g_.clear();
}

void
SingleStateMGXS::MakeSimple0(unsigned int num_groups, double sigma_t)
{
  Clear();

  num_groups_ = num_groups;
  sigma_t_.resize(num_groups, sigma_t);
  sigma_a_.resize(num_groups, sigma_t);

  ComputeDiffusionParameters();
}

void
SingleStateMGXS::MakeSimple1(unsigned int num_groups, double sigma_t, double c)
{
  Clear();

  num_groups_ = num_groups;
  sigma_t_.resize(num_groups, sigma_t);
  transfer_matrices_.emplace_back(num_groups, num_groups);

  // When multi-group, assign half the scattering cross section
  // to within-group scattering. The other half will be used for
  // up/down-scattering.
  auto& S = transfer_matrices_.back();
  double scale = (num_groups_ == 1) ? 1.0 : 0.5;
  S.SetDiagonal(std::vector<double>(num_groups, sigma_t * c * scale));

  // Set the up/down-scattering cross sections.
  // Summary:
  //     1) The half of groups with higher energies down-scatter to the next
  //        lowest energy group half the time and to the same group half the
  //        time.
  //     2) The half of groups with lower energies less the last group
  //        down-scatter to the next lowest energy group three quarters of the
  //        time and up-scatter to the next highest energy group one quarter
  //        of the time.
  //     3) The lowest energy group has the same form as 1).

  for (size_t g = 0; g < num_groups_; ++g)
  {
    // Down-scattering
    if (g > 0) S.Insert(g, g - 1, sigma_t * c * 0.5);

    // Up-scattering
    if (g > num_groups_ / 2)
    {
      if (g < num_groups_ - 1)
      {
        S.Insert(g, g - 1, sigma_t * c * 0.25);
        S.Insert(g, g + 1, sigma_t * c * 0.25);
      }
      else
        S.Insert(g, g - 1, sigma_t * c * 0.5);
    }
  } // for g

  ComputeAbsorption();
  ComputeDiffusionParameters();
}

void
SingleStateMGXS::MakeCombined(std::vector<std::pair<int, double>>& combinations)
{
  Clear();

  // Pickup all xs and make sure valid
  std::vector<std::shared_ptr<MultiGroupXS>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;
  double Nf_total = 0.0; // Total density of fissile materials

  // Loop over cross sections
  for (auto combo : combinations)
  {
    // Get the cross section from the stack
    std::shared_ptr<MultiGroupXS> xs;
    xs = Chi::GetStackItemPtr(Chi::multigroup_xs_stack, combo.first, std::string(__FUNCTION__));
    xsecs.push_back(xs);

    // Increment densities
    if (xs->IsFissionable())
    {
      is_fissionable_ = true;
      Nf_total += combo.second;
    }

    // Define and check number of groups
    if (xsecs.size() == 1) n_grps = xs->NumGroups();
    ChiLogicalErrorIf(xs->NumGroups() != n_grps,
                      "All cross sections being combined must have the same group structure.");

    // Increment number of precursors
    n_precs += xs->NumPrecursors();
  } // for cross section

  // Check that the fissile and precursor densities are greater than
  // machine precision. If this condition is not met, the material is assumed
  // to be either not fissile, have zero precursors, or both.
  if (Nf_total < 1.0e-12) is_fissionable_ = false;

  // Check to ensure that all fissionable cross sections contain either
  // prompt/delayed fission data or total fission data
  if (n_precs > 0)
    for (const auto& xs : xsecs)
      ChiLogicalErrorIf(xs->IsFissionable() and xs->NumPrecursors() == 0,
                        "If precursors are specified, all fissionable cross sections must "
                        "specify precursors.");

  // Initialize the data
  num_groups_ = n_grps;
  scattering_order_ = 0;
  num_precursors_ = n_precs;
  for (const auto& xs : xsecs)
    scattering_order_ = std::max(scattering_order_, xs->ScatteringOrder());

  // Init mandatory cross sections
  sigma_t_.assign(n_grps, 0.0);
  sigma_a_.assign(n_grps, 0.0);

  // Init transfer matrices only if at least one exists
  if (std::any_of(xsecs.begin(),
                  xsecs.end(),
                  [](const std::shared_ptr<MultiGroupXS>& x)
                  { return not x->TransferMatrices().empty(); }))
    transfer_matrices_.assign(scattering_order_ + 1, SparseMatrix(num_groups_, num_groups_));

  // Init fission data
  if (is_fissionable_)
  {
    sigma_f_.assign(n_grps, 0.0);
    nu_sigma_f_.assign(n_grps, 0.0);
    production_matrix_.assign(num_groups_, std::vector<double>(num_groups_, 0.0));

    // Init prompt/delayed fission data
    if (n_precs > 0)
    {
      nu_prompt_sigma_f_.assign(n_grps, 0.0);
      nu_delayed_sigma_f_.assign(n_grps, 0.0);
      precursors_.resize(n_precs);
    }
  }

  // Combine the data
  size_t precursor_count = 0;
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    // Atom density
    const auto N_i = combinations[x].second;

    // Fraction of fissile density
    double ff_i = 0.0;
    if (xsecs[x]->IsFissionable()) ff_i = N_i / Nf_total;

    // Combine cross sections
    const auto& sig_t = xsecs[x]->SigmaTotal();
    const auto& sig_a = xsecs[x]->SigmaAbsorption();
    const auto& sig_f = xsecs[x]->SigmaFission();
    const auto& nu_p_sig_f = xsecs[x]->NuPromptSigmaF();
    const auto& nu_d_sig_f = xsecs[x]->NuDelayedSigmaF();
    const auto& F = xsecs[x]->ProductionMatrix();

    // Here, raw cross sections are scaled by densities and spectra by
    // fractional densities. The latter is done to preserve a unit spectra.
    for (size_t g = 0; g < n_grps; ++g)
    {
      sigma_t_[g] += N_i * sig_t[g];
      sigma_a_[g] += N_i * sig_a[g];

      if (xsecs[x]->IsFissionable())
      {
        sigma_f_[g] += N_i * sig_f[g];
        nu_sigma_f_[g] += N_i * sig_f[g];
        for (size_t gp = 0; gp < num_groups_; ++gp)
          production_matrix_[g][gp] += N_i * F[g][gp];

        if (n_precs > 0)
        {
          nu_prompt_sigma_f_[g] += N_i * nu_p_sig_f[g];
          nu_delayed_sigma_f_[g] += N_i * nu_d_sig_f[g];
        }
      }
    } // for g

    // Combine precursor data
    // Here, all precursors across all materials are stored. The decay constants and delayed
    // spectrum are what they are, however, some special treatment must be given to the yields.
    // Because the yield dictates what fraction of delayed neutrons are produced from a
    // given family, the sum over all families must yield unity. To achieve this end, all
    // precursor yields must be scaled based on the fraction of the total density of materials
    // with precursors they make up.

    if (xsecs[x]->NumPrecursors() > 0)
    {
      const auto& precursors = xsecs[x]->Precursors();
      for (size_t j = 0; j < xsecs[x]->NumPrecursors(); ++j)
      {
        size_t count = precursor_count + j;
        const auto& precursor = precursors[j];
        precursors_[count].decay_constant = precursor.decay_constant;
        precursors_[count].fractional_yield = precursor.fractional_yield * ff_i;
        precursors_[count].emission_spectrum = precursor.emission_spectrum;
      } // for j

      precursor_count += xsecs[x]->NumPrecursors();
    }

    // Set inverse velocity data
    if (x == 0 and xsecs[x]->InverseVelocity().empty()) inv_velocity_ = xsecs[x]->InverseVelocity();
    ChiLogicalErrorIf(
      xsecs[x]->InverseVelocity() != inv_velocity_,
      "All cross sections being combined must have the same group-wise velocities.");

    // Combine transfer matrices

    // This step is somewhat tricky. The cross sections aren't guaranteed
    // to have the same sparsity patterns and therefore simply adding them
    // together has to take the sparse matrix's protection mechanisms into
    // account.

    if (not xsecs[x]->TransferMatrices().empty())
    {
      for (size_t m = 0; m < xsecs[x]->ScatteringOrder() + 1; ++m)
      {
        auto& Sm = transfer_matrices_[m];
        const auto& Sm_other = xsecs[x]->TransferMatrix(m);
        for (size_t g = 0; g < num_groups_; ++g)
        {
          const auto& cols = Sm_other.rowI_indices_[g];
          const auto& vals = Sm_other.rowI_values_[g];
          for (size_t t = 0; t < cols.size(); ++t)
            Sm.InsertAdd(g, t, vals[t] * N_i);
        }
      }
    }
  } // for cross sections

  ComputeDiffusionParameters();
}

void
SingleStateMGXS::MakeFromChiXSFile(const std::string& file_name)
{
  Clear();

  // Open Chi XS file
  std::ifstream file;
  file.open(file_name);
  ChiLogicalErrorIf(not file.is_open(), "Failed to open cross section file " + file_name + ".");
  log.Log() << "Reading Chi cross section file \"" << file_name << "\"\n";

  // Lambda for reading group structure data.
  auto ReadGroupStructure = [](const std::string& keyword,
                               std::vector<std::vector<double>>& destination,
                               const size_t n_grps,
                               std::ifstream& file,
                               std::istringstream& line_stream,
                               size_t& line_number)
  {
    destination.assign(n_grps, std::vector<double>(2, 0.0));

    std::string line;
    int group;
    double high, low;
    size_t count = 0;

    // Read the block
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    while (line != keyword + "_END")
    {
      // Get data from current line
      line_stream >> group >> high >> low;
      destination.at(group).at(0) = high;
      destination.at(group).at(1) = low;
      ChiLogicalErrorIf(count++ >= n_grps,
                        "Too many entries encountered when parsing group structure.\n"
                        "The expected number of entries is " +
                          std::to_string(n_grps) + ".");

      // Go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  // Lambda for reading vector data.
  auto Read1DData = [](const std::string& keyword,
                       std::vector<double>& destination,
                       const size_t n_entries,
                       std::ifstream& file,
                       std::istringstream& line_stream,
                       size_t& line_number)
  {
    destination.assign(n_entries, 0.0);

    std::string line;
    int i;
    double value;
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
      ChiLogicalErrorIf(count++ >= n_entries,
                        "Too many entries encountered when parsing 1D data.\n"
                        "The expected number of entries is " +
                          std::to_string(n_entries) + ".");

      // Go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  // Lambda for reading 2D data
  auto Read2DData = [](const std::string& keyword,
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
    double value;
    size_t i, j;

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
        if (entry_prefix == "G_PRECURSOR_VAL") destination.at(j).at(i) = value; // hack
        else
          destination.at(i).at(j) = value;
      }

      // Go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  // Lambda for reading transfer matrix data.
  auto ReadTransferMatrices = [](const std::string& keyword,
                                 std::vector<SparseMatrix>& destination,
                                 const size_t n_moms,
                                 const size_t n_grps,
                                 std::ifstream& file,
                                 std::istringstream& line_stream,
                                 size_t& line_number)
  {
    destination.assign(n_moms, SparseMatrix(n_grps, n_grps));

    std::string word, line;
    double value;
    size_t ell, group, gprime;

    // Read the block
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    while (line != keyword + "_END")
    {
      // Check that this line contains an entry
      line_stream >> word;
      if (word == "M_GPRIME_G_VAL")
      {
        // Get data from current line
        line_stream >> ell >> gprime >> group >> value;
        destination.at(ell).Insert(group, gprime, value);
      }

      // Go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  // Lambda for checking for all non-negative values.
  auto IsNonNegative = [](const std::vector<double>& vec)
  { return not std::any_of(vec.begin(), vec.end(), [](double x) { return x < 0.0; }); };

  /// Lambda for checking for all strictly positive values.
  auto IsPositive = [](const std::vector<double>& vec)
  { return not std::any_of(vec.begin(), vec.end(), [](double x) { return x <= 0.0; }); };

  /// Lambda for checking for any non-zero values.
  auto HasNonZero = [](const std::vector<double>& vec)
  { return std::any_of(vec.begin(), vec.end(), [](double x) { return x > 0.0; }); };

  // Read the Chi XS file
  // TODO: Determine whether or not to allow specification of a
  //       data block without any data. Currently, if a data block
  //       is specified and no values are present, the std::any_of
  //       checks will evaluate false if expected data is not present.

  std::vector<double> decay_constants;
  std::vector<double> fractional_yields;
  std::vector<std::vector<double>> emission_spectra;
  std::vector<double> nu, nu_prompt, nu_delayed, beta;
  std::vector<double> chi, chi_prompt;

  std::string word, line;
  size_t line_number = 0;
  while (std::getline(file, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    // Parse number of groups
    if (word == "NUM_GROUPS")
    {
      int n_groups;
      line_stream >> n_groups;
      ChiLogicalErrorIf(n_groups <= 0, "The number of energy groups must be positive.");
      num_groups_ = n_groups;
    }

    // Parse the number of scattering moments
    if (word == "NUM_MOMENTS")
    {
      int n_moments;
      line_stream >> n_moments;
      ChiLogicalErrorIf(n_moments < 0, "The number of scattering moments must be non-negative.");
      scattering_order_ = std::max(0, n_moments - 1);
    }

    // Parse the number of precursors species
    if (word == "NUM_PRECURSORS")
    {
      int n_prec;
      line_stream >> n_prec;
      ChiLogicalErrorIf(n_prec < 0,
                        "The number of delayed neutron precursors must be non-negative.");
      num_precursors_ = n_prec;
      precursors_.resize(num_precursors_);
    }

    // Parse nuclear data
    try
    {
      auto& ln = line_number;
      auto& ls = line_stream;
      auto& f = file;
      auto& fw = word;

      //
      // Group Structure Data
      //

      if (fw == "GROUP_STRUCTURE_BEGIN")
        ReadGroupStructure("GROUP_STRUCTURE", e_bounds_, num_groups_, f, ls, ln);
      if (fw == "INV_VELOCITY_BEGIN")
      {
        Read1DData("INV_VELOCITY", inv_velocity_, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not IsPositive(inv_velocity_),
                          "Only positive inverse velocity values are permitted.");
      }
      if (fw == "VELOCITY_BEGIN" and inv_velocity_.empty())
      {
        Read1DData("VELOCITY", inv_velocity_, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not IsPositive(inv_velocity_),
                          "Only positive velocity values are permitted.");

        // Compute inverse velocity
        for (size_t g = 0; g < num_groups_; ++g)
          inv_velocity_[g] = 1.0 / inv_velocity_[g];
      }

      //
      // Read cross section data
      //

      if (fw == "SIGMA_T_BEGIN")
      {
        Read1DData("SIGMA_T", sigma_t_, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not IsNonNegative(sigma_t_),
                          "Only non-negative total cross section values are permitted.");
      } // if sigma_t

      if (fw == "SIGMA_A_BEGIN")
      {
        Read1DData("SIGMA_A", sigma_a_, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not IsNonNegative(sigma_a_),
                          "Only non-negative absorption cross section values are permitted.");
      } // if sigma_a

      if (fw == "SIGMA_F_BEGIN")
      {
        Read1DData("SIGMA_F", sigma_f_, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not IsNonNegative(sigma_f_),
                          "Only non-negative fission cross section values are permitted.");
        if (not HasNonZero(sigma_f_))
        {
          log.Log0Warning() << "The fission cross section specified in "
                            << "\"" << file_name << "\" is uniformly zero... Clearing it.";
          sigma_f_.clear();
        }
      } // if sigma_f

      if (fw == "NU_SIGMA_F_BEGIN")
      {
        Read1DData("NU_SIGMA_F", nu_sigma_f_, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(
          not IsNonNegative(nu_sigma_f_),
          "Only non-negative total fission multiplication cross section values are permitted.");
        if (not HasNonZero(nu_sigma_f_))
        {
          log.Log0Warning() << "The production cross-section specified in "
                            << "\"" << file_name << "\" is uniformly zero... Clearing it.";
          nu_sigma_f_.clear();
        }
      } // if nu_sigma_f

      //
      // Read neutrons per fission data
      //

      if (fw == "NU_BEGIN")
      {
        Read1DData("NU", nu, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(
          not std::all_of(nu.begin(), nu.end(), [](double x) { return x == 0.0 or x > 1.0; }),
          "Total fission neutron yield values must be either zero, or greater than one.");
        if (not HasNonZero(nu))
        {
          log.Log0Warning() << "The total fission neutron yield specified in "
                            << "\"" << file_name << "\" is uniformly zero... Clearing it.";
          nu.clear();
        }

        // Compute prompt/delayed nu, if needed
        if (num_precursors_ > 0 and not nu.empty() and not beta.empty() and nu_prompt.empty() and
            nu_delayed.empty())
        {
          nu_prompt.assign(num_groups_, 0.0);
          nu_delayed.assign(num_groups_, 0.0);
          for (size_t g = 0; g < num_groups_; ++g)
          {
            nu_prompt[g] = (1.0 - beta[g]) * nu[g];
            nu_delayed[g] = beta[g] * nu[g];
          }
        }
      } // if nu

      if (fw == "NU_PROMPT_BEGIN")
      {
        Read1DData("NU_PROMPT", nu_prompt, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not std::all_of(nu_prompt.begin(),
                                          nu_prompt.end(),
                                          [](double x) { return x == 0.0 or x > 1.0; }),
                          "Average prompt fission neutron yield values must be either zero, "
                          "or greater than one.");
        if (not HasNonZero(nu_prompt))
        {
          log.Log0Warning() << "The prompt fission neutron yield specified in "
                            << "\"" << file_name << "\" is uniformly zero... Clearing it.";
          nu_prompt.clear();
        }
      } // if nu_prompt

      if (fw == "NU_DELAYED_BEGIN")
      {
        Read1DData("NU_DELAYED", nu_delayed, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not IsNonNegative(nu_delayed),
                          "Average delayed fission neutron yield values "
                          "must be non-negative.");
        if (not HasNonZero(nu_delayed))
        {
          log.Log0Warning() << "The delayed fission neutron yield specified in "
                            << "\"" << file_name << "\" is uniformly zero... Clearing it.";
          nu_prompt.clear();
        }
      } // if nu_delayed

      if (fw == "BETA_BEGIN")
      {
        Read1DData("BETA", beta, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(
          not std::all_of(beta.begin(), beta.end(), [](double x) { return x >= 0.0 and x <= 1.0; }),
          "Delayed neutron fraction values must be in the range [0.0, 1.0].");
        if (not HasNonZero(beta))
        {
          log.Log0Warning() << "The delayed neutron fraction specified in "
                            << "\"" << file_name << "\" is uniformly zero... Clearing it.";
          beta.clear();
        }

        // Compute prompt/delayed nu, if needed
        if (num_precursors_ > 0 and not nu.empty() and not beta.empty() and nu_prompt.empty() and
            nu_delayed.empty())
        {
          nu_prompt.assign(num_groups_, 0.0);
          nu_delayed.assign(num_groups_, 0.0);
          for (unsigned int g = 0; g < num_groups_; ++g)
          {
            nu_prompt[g] = (1.0 - beta[g]) * nu[g];
            nu_delayed[g] = beta[g] * nu[g];
          }
        }
      } // if beta

      //
      // Read fission/emission spectra
      //

      if (fw == "CHI_BEGIN")
      {
        Read1DData("CHI", chi, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(
          not HasNonZero(chi),
          "The steady-state fission spectrum must have at least one non-zero value.");
        ChiLogicalErrorIf(not IsNonNegative(chi),
                          "The steady-state fission spectrum must be non-negative.");

        // Normalizing
        const auto sum = std::accumulate(chi.begin(), chi.end(), 0.0);
        std::transform(chi.begin(), chi.end(), chi.begin(), [sum](double& x) { return x / sum; });
      } // if chi

      if (fw == "CHI_PROMPT_BEGIN")
      {
        Read1DData("CHI_PROMPT", chi_prompt, num_groups_, f, ls, ln);
        ChiLogicalErrorIf(not HasNonZero(chi_prompt),
                          "The prompt fission spectrum must have at least one non-zero value.");
        ChiLogicalErrorIf(not IsNonNegative(chi_prompt),
                          "The prompt fission spectrum must be non-negative.");

        // Normalizing
        const auto sum = std::accumulate(chi_prompt.begin(), chi_prompt.end(), 0.0);
        std::transform(chi_prompt.begin(),
                       chi_prompt.end(),
                       chi_prompt.begin(),
                       [sum](double& x) { return x / sum; });

      } // if prompt chi

      if (num_precursors_ > 0 and fw == "CHI_DELAYED_BEGIN")
      {
        // TODO: Should the be flipped to PRECURSOR_G_VAL?
        Read2DData("CHI_DELAYED",
                   "G_PRECURSOR_VAL",
                   emission_spectra,
                   num_precursors_,
                   num_groups_,
                   f,
                   ls,
                   ln);

        for (size_t j = 0; j < num_precursors_; ++j)
        {
          ChiLogicalErrorIf(not HasNonZero(emission_spectra[j]),
                            "Delayed emission spectrum for precursor " + std::to_string(j) +
                              " must have at least one non-zero value.");
          ChiLogicalErrorIf(not IsNonNegative(emission_spectra[j]),
                            "Delayed emission spectrum for precursor " + std::to_string(j) +
                              " must be non-negative.");

          // normalizing
          const auto sum =
            std::accumulate(emission_spectra[j].begin(), emission_spectra[j].end(), 0.0);
          std::transform(emission_spectra[j].begin(),
                         emission_spectra[j].end(),
                         emission_spectra[j].begin(),
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
          Read1DData("PRECURSOR_DECAY_CONSTANTS", decay_constants, num_precursors_, f, ls, ln);
          ChiLogicalErrorIf(not IsPositive(decay_constants),
                            "Delayed neutron precursor decay constants must be positive.");
        } // if decay constants

        if (fw == "PRECURSOR_FRACTIONAL_YIELDS_BEGIN")
        {
          Read1DData("PRECURSOR_FRACTIONAL_YIELDS", fractional_yields, num_precursors_, f, ls, ln);
          ChiLogicalErrorIf(
            not HasNonZero(fractional_yields),
            "Delayed neutron precursor fractional yields must contain at least one non-zero.");
          ChiLogicalErrorIf(
            not std::all_of(fractional_yields.begin(),
                            fractional_yields.end(),
                            [](double x) { return x >= 0.0 and x <= 1.0; }),
            "Delayed neutron precursor fractional yields must be in the range [0.0, 1.0].");

          // Normalizing
          const auto sum = std::accumulate(fractional_yields.begin(), fractional_yields.end(), 0.0);
          std::transform(fractional_yields.begin(),
                         fractional_yields.end(),
                         fractional_yields.begin(),
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
      throw std::runtime_error("Error reading Chi cross section file "
                               "\"" +
                               file_name + "\".\n" + "Line number " + std::to_string(line_number) +
                               "\n" + err.what());
    }

    catch (const std::logic_error& err)
    {
      throw std::logic_error("Error reading Chi cross section file "
                             "\"" +
                             file_name + "\".\n" + "Line number " + std::to_string(line_number) +
                             "\n" + err.what());
    }

    catch (...)
    {
      throw std::runtime_error("Unknown error encountered.");
    }

    word = "";
  } // while not EOF, read each lines
  file.close();

  if (sigma_a_.empty()) ComputeAbsorption();
  ComputeDiffusionParameters();

  //
  // Compute and check fission data
  //

  // Determine if the material is fissionable
  is_fissionable_ =
    not sigma_f_.empty() or not nu_sigma_f_.empty() or not production_matrix_.empty();

  // Check and set the fission data
  if (is_fissionable_)
  {
    // Check vector data inputs
    if (production_matrix_.empty())
    {
      // Check for non-delayed fission neutron yield data
      ChiLogicalErrorIf(nu.empty() and nu_prompt.empty(),
                        "Either the total or prompt fission neutron yield must be specified "
                        "for fissionable materials.");
      ChiLogicalErrorIf(not nu.empty() and not nu_prompt.empty(),
                        "Ambiguous fission neutron yield. Only one of the total and prompt "
                        "fission neutron yield should be specified.");

      // Check for fission spectrum data
      ChiLogicalErrorIf(chi.empty() and chi_prompt.empty(),
                        "Either the steady-state or prompt fission spectrum must be specified "
                        "for fissionable materials.");
      ChiLogicalErrorIf(not chi.empty() and not chi_prompt.empty(),
                        "Ambiguous fission spectrum data. Only one of the steady-state and "
                        "prompt fission spectrum should be specified.");

      // Check for compatibility
      if ((not nu.empty() and chi.empty()) or (nu.empty() and not chi.empty()) or
          (not nu_prompt.empty() and chi_prompt.empty()) or
          (nu_prompt.empty() and not chi_prompt.empty()))
        ChiLogicalError("Ambiguous fission data. Either the total fission neutron yield with the "
                        "steady-state fission spectrum or the prompt fission neutron yield with "
                        "the prompt fission spectrum should be specified.");

      // Initialize total fission neutron yield from prompt
      if (not nu_prompt.empty()) nu = nu_prompt;

      // Check delayed neutron data
      if (num_precursors_ > 0)
      {
        // Check that decay data was specified
        ChiLogicalErrorIf(decay_constants.empty(),
                          "Precursor decay constants are required when precursors are specified.");

        // Check that yield data was specified
        ChiLogicalErrorIf(fractional_yields.empty(),
                          "Precursor yields are required when precursors are specified.");

        // Check that prompt data was specified
        ChiLogicalErrorIf(chi_prompt.empty() or nu_prompt.empty(),
                          "Both the prompt fission spectrum and prompt fission neutron yield "
                          "must be specified when delayed neutron precursors are specified.");

        // Check that delayed neutron production and emission spectra were specified
        ChiLogicalErrorIf(nu_delayed.empty() or
                            std::any_of(emission_spectra.begin(),
                                        emission_spectra.end(),
                                        [](const std::vector<double>& x) { return x.empty(); }),
                          "Both the delay emission spectra and delayed fission neutron yield "
                          "must be specified when precursors are specified.");

        // Add delayed fission neutron yield to total
        for (size_t g = 0; g < num_groups_; ++g)
          nu[g] += nu_delayed[g];

        // Add data to precursor structs
        for (size_t j = 0; j < num_precursors_; ++j)
        {
          precursors_[j].decay_constant = decay_constants[j];
          precursors_[j].fractional_yield = fractional_yields[j];
          precursors_[j].emission_spectrum = emission_spectra[j];
        }
      }

      // Compute fission cross section
      if (sigma_f_.empty() and not nu_sigma_f_.empty())
      {
        sigma_f_ = nu_sigma_f_;
        for (size_t g = 0; g < num_groups_; ++g)
          if (nu_sigma_f_[g] > 0.0) sigma_f_[g] /= nu[g];
      }

      // Compute total production cross section
      nu_sigma_f_ = sigma_f_;
      for (size_t g = 0; g < num_groups_; ++g)
        nu_sigma_f_[g] *= nu[g];

      // Compute prompt production cross section
      if (not nu_prompt.empty())
      {
        nu_prompt_sigma_f_ = sigma_f_;
        for (size_t g = 0; g < num_groups_; ++g)
          nu_prompt_sigma_f_[g] *= nu_prompt[g];
      }

      // Compute delayed production cross section
      if (not nu_delayed.empty())
      {
        nu_delayed_sigma_f_ = sigma_f_;
        for (size_t g = 0; g < num_groups_; ++g)
          nu_delayed_sigma_f_[g] *= nu_delayed[g];
      }

      // Compute production matrix
      const auto fis_spec = not chi_prompt.empty() ? chi_prompt : chi;
      const auto nu_sigma_f = not nu_prompt.empty() ? nu_prompt_sigma_f_ : nu_sigma_f_;

      production_matrix_.resize(num_groups_);
      for (size_t g = 0; g < num_groups_; ++g)
        for (size_t gp = 0.0; gp < num_groups_; ++gp)
          production_matrix_[g].push_back(fis_spec[g] * nu_sigma_f[gp]);
    } // if production_matrix empty

    else
    {
      // TODO: Develop an implementation for multi-particle delayed neutron data.
      //       The primary challenge in this is that different precursor species exist for
      //       neutron-induced fission than for photo-fission.

      ChiLogicalErrorIf(num_precursors_ > 0,
                        "Currently, production matrix specification is not allowed when "
                        "delayed neutrons are present.");

      // Check for fission cross sections
      ChiLogicalErrorIf(sigma_f_.empty(),
                        "When a production matrix is specified, it must "
                        "be accompanied with a fission cross section.");

      // Compute production cross section
      nu_sigma_f_.assign(num_groups_, 0.0);
      for (size_t g = 0; g < num_groups_; ++g)
        for (size_t gp = 0; gp < num_groups_; ++gp)
          nu_sigma_f_[gp] += production_matrix_[g][gp];

      // Check for reasonable fission neutron yield
      nu = nu_sigma_f_;
      for (size_t g = 0; g < num_groups_; ++g)
        if (sigma_f_[g] > 0.0) nu[g] /= sigma_f_[g];

      ChiLogicalErrorIf(
        not IsNonNegative(nu),
        "The production matrix implies an invalid negative average fission neutron yield.");
      ChiLogicalErrorIf(
        not std::all_of(nu.begin(), nu.end(), [](double x) { return x == 0.0 and x > 1.0; }),
        "Incompatible fission data encountered. The computed nu is not either zero or "
        "greater than one.");

      if (std::any_of(nu.begin(), nu.end(), [](double x) { return x > 8.0; }))
        log.Log0Warning() << "A computed nu of greater than 8.0 was encountered. ";
    }

    ChiLogicalErrorIf(
      sigma_f_.empty(),
      "Fissionable materials are required to have a defined fission cross section.");
  } // if fissionable

  // Clear fission data if not fissionable
  else
  {
    sigma_f_.clear();
    nu_sigma_f_.clear();
    nu_prompt_sigma_f_.clear();
    nu_delayed_sigma_f_.clear();
    production_matrix_.clear();
    precursors_.clear();
  } // if not fissionable
}

void
SingleStateMGXS::ComputeAbsorption()
{
  sigma_a_.assign(num_groups_, 0.0);

  // Compute for a pure absorber
  if (transfer_matrices_.empty())
    for (size_t g = 0; g < num_groups_; ++g)
      sigma_a_[g] = sigma_t_[g];

  // Estimate from a transfer matrix
  else
  {
    log.Log0Warning() << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices_.front();
    for (size_t g = 0; g < num_groups_; ++g)
    {
      // Estimate the scattering cross section
      double sigma_s = 0.0;
      for (size_t row = 0; row < S0.NumRows(); ++row)
      {
        const auto& cols = S0.rowI_indices_[row];
        const auto& vals = S0.rowI_values_[row];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sigma_s += vals[t];
            break;
          }
      }
      sigma_a_[g] = sigma_t_[g] - sigma_s;

      // TODO: Should negative absorption be allowed?
      if (sigma_a_[g] < 0.0)
        log.Log0Warning() << "Negative absorption cross section encountered in group " << g
                          << " when estimating from the transfer matrices";
    } // for g
  }   // if scattering present
}

void
SingleStateMGXS::ComputeDiffusionParameters()
{
  if (diffusion_initialized_) return;

  // Initialize diffusion data
  sigma_tr_.resize(num_groups_, 0.0);
  diffusion_coeff_.resize(num_groups_, 1.0);
  sigma_s_gtog_.resize(num_groups_, 0.0);
  sigma_r_.resize(num_groups_, 0.0);

  // Perform computations group-wise
  const auto& S = transfer_matrices_;
  for (size_t g = 0; g < num_groups_; ++g)
  {
    // Determine transport correction
    double sigma_1 = 0.0;
    if (S.size() > 1)
    {
      for (size_t gp = 0; gp < num_groups_; ++gp)
      {
        const auto& cols = S[1].rowI_indices_[gp];
        const auto& vals = S[1].rowI_values_[gp];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sigma_1 += vals[t];
            break;
          }
      } // for gp
    }   // if moment 1 available

    // Compute transport cross section
    if (sigma_1 >= sigma_t_[g])
    {
      log.Log0Warning() << "Transport corrected diffusion coefficient failed for group " << g
                        << " in call to " << __FUNCTION__ << ". "
                        << "sigma_t=" << sigma_t_[g] << " sigma_1=" << sigma_1 << ". "
                        << "Setting sigma_1 to zero for group " << g << ".";
      sigma_1 = 0.0;
    }
    sigma_tr_[g] = sigma_t_[g] - sigma_1;

    // Compute the diffusion coefficient
    // Cap the value for when sigma_t - sigma_1 is near zero
    diffusion_coeff_[g] = std::fmin(1.0e12, 1.0 / 3.0 / (sigma_t_[g] - sigma_1));

    // Determine within group scattering
    if (not S.empty())
    {
      const auto& cols = S[0].rowI_indices_[g];
      const auto& vals = S[0].rowI_values_[g];
      for (size_t t = 0; t < cols.size(); ++t)
        if (cols[t] == g)
        {
          sigma_s_gtog_[g] = vals[t];
          break;
        }
    }

    // Compute removal cross section
    sigma_r_[g] = std::max(0.0, sigma_t_[g] - sigma_s_gtog_[g]);
  } // for g

  diffusion_initialized_ = true;
}

} // namespace opensn
