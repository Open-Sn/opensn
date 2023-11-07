#include "framework/physics/PhysicsMaterial/MultiGroupXS/single_state_mgxs.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include <algorithm>
#include <string>
#include <numeric>

namespace chi_physics
{

void
SingleStateMGXS::Clear()
{
  num_groups_ = 0;
  scattering_order_ = 0;
  num_precursors_ = 0;
  is_fissionable_ = false;

  sigma_t_.clear();
  sigma_f_.clear();
  sigma_a_.clear();

  nu_sigma_f_.clear();
  nu_prompt_sigma_f_.clear();
  nu_delayed_sigma_f_.clear();

  inv_velocity_.clear();

  transfer_matrices_.clear();
  production_matrix_.clear();

  precursors_.clear();

  // Diffusion quantities
  diffusion_initialized_ = false;
  diffusion_coeff_.clear();
  sigma_removal_.clear();
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

  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    // downscattering
    if (g > 0) S.Insert(g, g - 1, sigma_t * c * 0.5);

    // upscattering
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

  // pickup all xs and make sure valid
  std::vector<std::shared_ptr<MultiGroupXS>> xsecs;
  xsecs.reserve(combinations.size());

  unsigned int n_grps = 0;
  unsigned int n_precs = 0;
  double Nf_total = 0.0; // total density of fissile materials

  // loop over cross sections
  for (auto combo : combinations)
  {
    // get the cross section from the lua stack
    std::shared_ptr<MultiGroupXS> xs;
    xs = Chi::GetStackItemPtr(Chi::multigroup_xs_stack, combo.first, std::string(__FUNCTION__));
    xsecs.push_back(xs);

    // increment densities
    if (xs->IsFissionable())
    {
      is_fissionable_ = true;
      Nf_total += combo.second;
    }

    // define and check number of groups
    if (xsecs.size() == 1) n_grps = xs->NumGroups();
    else if (xs->NumGroups() != n_grps)
      throw std::logic_error("Incompatible cross sections encountered.\n"
                             "All cross sections being combined must have the "
                             "same number of energy groups.");

    // increment number of precursors
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
      if (xs->IsFissionable() and xs->NumPrecursors() == 0)
        throw std::logic_error("Incompatible cross sections encountered.\n"
                               "If any fissionable cross sections specify delayed neutron "
                               "data, all fissionable cross sections must specify delayed "
                               "neutron data.");

  // Initialize the data
  num_groups_ = n_grps;
  num_precursors_ = n_precs;
  scattering_order_ = 0;
  for (const auto& xs : xsecs)
    scattering_order_ = std::max(scattering_order_, xs->ScatteringOrder());

  // mandatory cross sections
  sigma_t_.assign(n_grps, 0.0);
  sigma_a_.assign(n_grps, 0.0);

  // init transfer matrices only if at least one exists
  using XSPtr = MultiGroupXSPtr;
  if (std::any_of(xsecs.begin(),
                  xsecs.end(),
                  [](const XSPtr& x) { return not x->TransferMatrices().empty(); }))
    transfer_matrices_.assign(scattering_order_ + 1,
                              chi_math::SparseMatrix(num_groups_, num_groups_));

  // init fission data
  if (is_fissionable_)
  {
    sigma_f_.assign(n_grps, 0.0);
    nu_sigma_f_.assign(n_grps, 0.0);
    production_matrix_.assign(num_groups_, std::vector<double>(num_groups_, 0.0));

    // init prompt/delayed fission data
    if (n_precs > 0)
    {
      nu_prompt_sigma_f_.assign(n_grps, 0.0);
      nu_delayed_sigma_f_.assign(n_grps, 0.0);
      precursors_.resize(n_precs);
    }
  }

  // Combine the data
  unsigned int precursor_count = 0;
  for (size_t x = 0; x < xsecs.size(); ++x)
  {
    // atom density
    double N_i = combinations[x].second;

    // fraction of fissile density
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
    for (unsigned int g = 0; g < n_grps; ++g)
    {
      sigma_t_[g] += sig_t[g] * N_i;
      sigma_a_[g] += sig_a[g] * N_i;

      if (xsecs[x]->IsFissionable())
      {
        sigma_f_[g] += sig_f[g] * N_i;
        nu_sigma_f_[g] += sig_f[g] * N_i;
        for (unsigned int gp = 0; gp < num_groups_; ++gp)
          production_matrix_[g][gp] += F[g][gp] * N_i;

        if (n_precs > 0)
        {
          nu_prompt_sigma_f_[g] += nu_p_sig_f[g] * N_i;
          nu_delayed_sigma_f_[g] += nu_d_sig_f[g] * N_i;
        }
      }
    } // for g

    // Combine precursor data

    // Here, all precursors across all materials are stored. The decay
    // constants and delayed spectrum are what they are, however, some
    // special treatment must be given to the yields. Because the yield
    // tells us what fraction of delayed neutrons are produced from a
    // given family, the sum over all families must yield unity. To
    // achieve this end, we must scale all precursor yields based on
    // the fraction of the total density of materials with precursors
    // they make up.

    if (xsecs[x]->NumPrecursors() > 0)
    {
      const auto& precursors = xsecs[x]->Precursors();
      for (unsigned int j = 0; j < xsecs[x]->NumPrecursors(); ++j)
      {
        unsigned int count = precursor_count + j;
        const auto& precursor = precursors[j];
        precursors_[count].decay_constant = precursor.decay_constant;
        precursors_[count].fractional_yield = precursor.fractional_yield * ff_i;
        precursors_[count].emission_spectrum = precursor.emission_spectrum;
      } // for j

      precursor_count += xsecs[x]->NumPrecursors();
    }

    // Set inverse velocity data
    if (x == 0 && !xsecs[x]->InverseVelocity().empty()) inv_velocity_ = xsecs[x]->InverseVelocity();
    else if (xsecs[x]->InverseVelocity() != inv_velocity_)
      throw std::logic_error("Invalid cross sections encountered.\n"
                             "All cross sections being combined must share a group "
                             "structure. This implies that the inverse speeds for "
                             "each of the cross sections must be equivalent.");

    // Combine transfer matrices

    // This step is somewhat tricky. The cross sections aren't guaranteed
    // to have the same sparsity patterns and therefore simply adding them
    // together has to take the sparse matrix's protection mechanisms into
    // account.

    if (not xsecs[x]->TransferMatrices().empty())
    {
      for (unsigned int m = 0; m < xsecs[x]->ScatteringOrder() + 1; ++m)
      {
        auto& Sm = transfer_matrices_[m];
        const auto& Sm_other = xsecs[x]->TransferMatrix(m);
        for (unsigned int g = 0; g < num_groups_; ++g)
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

  Chi::log.Log() << "Reading Chi cross section file \"" << file_name << "\"\n";

  // opens and checks if open
  std::ifstream file;
  file.open(file_name);
  if (not file.is_open())
    throw std::runtime_error("Failed to open Chi cross section file "
                             "\"" +
                             file_name + "\" in call to " + std::string(__FUNCTION__));

  // Define utility functions for parsing

  /// Lambda for reading group structure data.
  auto ReadGroupStructure = [](const std::string& keyword,
                               std::vector<std::vector<double>>& destination,
                               const unsigned int G, // # of groups
                               std::ifstream& file,
                               std::istringstream& line_stream,
                               unsigned int& line_number)
  {
    // init storage
    destination.assign(G, std::vector<double>(2, 0.0));

    // bookkeeping
    std::string line;
    int group;
    double high;
    double low;
    unsigned int count = 0;

    // read the block
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    while (line != keyword + "_END")
    {
      // get data from current line
      line_stream >> group >> high >> low;
      destination.at(group).at(0) = high;
      destination.at(group).at(1) = low;
      if (count++ >= G)
        throw std::runtime_error("Too many entries encountered when parsing "
                                 "group structure.\nThe expected number of entries "
                                 "is " +
                                 std::to_string(G) + ".");

      // go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  /// Lambda for reading vector data.
  auto Read1DData = [](const std::string& keyword,
                       std::vector<double>& destination,
                       const unsigned int N,
                       std::ifstream& file,
                       std::istringstream& line_stream,
                       unsigned int& line_number)
  {
    // init storage
    destination.assign(N, 0.0);

    // bookkeeping
    std::string line;
    int i;
    double value;
    unsigned int count = 0;

    // read the bloc
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    while (line != keyword + "_END")
    {
      // get data from current line
      line_stream >> i >> value;
      destination.at(i) = value;
      if (count++ >= N)
        throw std::runtime_error("To many entries encountered when parsing "
                                 "1D data.\nThe expected number of entries is " +
                                 std::to_string(N));

      // go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  /// Lambda for reading 2D data
  auto Read2DData = [](const std::string& keyword,
                       const std::string& entry_prefix,
                       std::vector<std::vector<double>>& destination,
                       const unsigned int n_rows,
                       const unsigned int n_cols,
                       std::ifstream& file,
                       std::istringstream& line_stream,
                       unsigned int& line_number)
  {
    // init storage
    destination.clear();
    for (unsigned int i = 0; i < n_rows; ++i)
      destination.emplace_back(n_cols, 0.0);

    // bookkeeping
    std::string word, line;
    double value;
    unsigned int i;
    unsigned int j;

    // read the block
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    while (line != keyword + "_END")
    {
      // check that this line contains an entry
      line_stream >> word;
      if (word == entry_prefix)
      {
        // get data from current line
        line_stream >> i >> j >> value;
        if (entry_prefix == "G_PRECURSOR_VAL") destination.at(j).at(i) = value; // hack
        else
          destination.at(i).at(j) = value;
      }

      // go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  /// Lambda for reading transfer matrix data.
  auto ReadTransferMatrices = [](const std::string& keyword,
                                 std::vector<chi_math::SparseMatrix>& destination,
                                 const unsigned int M, // # of moments
                                 const unsigned int G, // # of groups
                                 std::ifstream& file,
                                 std::istringstream& line_stream,
                                 unsigned int& line_number)
  {
    // init storage
    destination.clear();
    for (unsigned int i = 0; i < M; ++i)
      destination.emplace_back(G, G);

    // bookkeeping
    std::string word, line;
    double value;
    unsigned int ell;
    unsigned int group;
    unsigned int gprime;

    // read the block
    std::getline(file, line);
    line_stream = std::istringstream(line);
    ++line_number;
    while (line != keyword + "_END")
    {
      // check that this line contains an entry
      line_stream >> word;
      if (word == "M_GPRIME_G_VAL")
      {
        // get data from current line
        line_stream >> ell >> gprime >> group >> value;
        destination.at(ell).Insert(group, gprime, value);
      }

      // go to next line
      std::getline(file, line);
      line_stream = std::istringstream(line);
      ++line_number;
    }
  };

  /// Lambda for checking for all non-negative values.
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
  unsigned int line_number = 0;
  while (std::getline(file, line))
  {
    std::istringstream line_stream(line);
    line_stream >> word;

    // parse number of groups
    if (word == "NUM_GROUPS")
    {
      int G;
      line_stream >> G;
      if (G <= 0)
        throw std::logic_error("The specified number of energy groups "
                               "must be strictly positive.");
      num_groups_ = G;
    }

    // parse the number of scattering moments
    if (word == "NUM_MOMENTS")
    {
      int M;
      line_stream >> M;
      if (M < 0)
        throw std::logic_error("The specified number of scattering moments "
                               "must be non-negative.");
      scattering_order_ = std::max(0, M - 1);
    }

    // parse the number of precursors species
    if (word == "NUM_PRECURSORS")
    {
      int J;
      line_stream >> J;
      if (J < 0)
        throw std::logic_error("The specified number of delayed neutron "
                               "precursors must be non-negative.");
      num_precursors_ = J;
      precursors_.resize(num_precursors_);
    }

    // parse nuclear data
    try
    {
      auto& ln = line_number;
      auto& ls = line_stream;
      auto& f = file;
      auto& fw = word;

      // Group Structure Data

      if (fw == "GROUP_STRUCTURE_BEGIN")
        ReadGroupStructure("GROUP_STRUCTURE", e_bounds_, num_groups_, f, ls, ln);

      if (fw == "INV_VELOCITY_BEGIN")
      {
        Read1DData("INV_VELOCITY", inv_velocity_, num_groups_, f, ls, ln);

        if (not IsPositive(inv_velocity_))
          throw std::logic_error("Invalid inverse velocity value encountered.\n"
                                 "Only strictly positive values are permitted.");
      }
      if (fw == "VELOCITY_BEGIN" and inv_velocity_.empty())
      {
        Read1DData("VELOCITY", inv_velocity_, num_groups_, f, ls, ln);

        if (not IsPositive(inv_velocity_))
          throw std::logic_error("Invalid velocity value encountered.\n"
                                 "Only strictly positive values are permitted.");

        // compute inverse
        for (unsigned int g = 0; g < num_groups_; ++g)
          inv_velocity_[g] = 1.0 / inv_velocity_[g];
      }

      // Cross Section Data

      if (fw == "SIGMA_T_BEGIN")
      {
        Read1DData("SIGMA_T", sigma_t_, num_groups_, f, ls, ln);

        if (not IsNonNegative(sigma_t_))
          throw std::logic_error("Invalid total cross section value encountered.\n"
                                 "Negative values are not permitted.");
      } // if sigma_t

      if (fw == "SIGMA_A_BEGIN")
      {
        Read1DData("SIGMA_A", sigma_a_, num_groups_, f, ls, ln);

        if (not IsNonNegative(sigma_a_))
          throw std::logic_error("Invalid absorption cross section value encountered.\n"
                                 "Negative values are not permitted.");
      } // if sigma_a

      if (fw == "SIGMA_F_BEGIN")
      {
        Read1DData("SIGMA_F", sigma_f_, num_groups_, f, ls, ln);

        if (not HasNonZero(sigma_f_))
        {
          Chi::log.Log0Warning() << "The fission cross section specified in "
                                 << "\"" << file_name << "\" is uniformly zero..."
                                 << "Clearing it.";
          sigma_f_.clear();
        }
        else if (not IsNonNegative(sigma_f_))
          throw std::logic_error("Invalid fission cross section value encountered.\n"
                                 "Negative values are not permitted.");
      } // if sigma_f

      if (fw == "NU_SIGMA_F_BEGIN")
      {
        Read1DData("NU_SIGMA_F", nu_sigma_f_, num_groups_, f, ls, ln);

        if (not HasNonZero(nu_sigma_f_))
        {
          Chi::log.Log0Warning() << "The production cross-section specified in "
                                 << "\"" << file_name << "\" is uniformly zero..."
                                 << "Clearing it.";
          nu_sigma_f_.clear();
        }
        else if (not IsNonNegative(nu_sigma_f_))
          throw std::logic_error("Invalid production cross section value encountered.\n"
                                 "Negative values are not permitted.");
      } // if nu_sigma_f

      // Neutrons Per Fission

      if (fw == "NU_BEGIN")
      {
        Read1DData("NU", nu, num_groups_, f, ls, ln);

        if (not HasNonZero(nu))
        {
          Chi::log.Log0Warning() << "The total fission neutron yield specified in "
                                 << "\"" << file_name << "\" is uniformly zero..."
                                 << "Clearing it.";
          nu.clear();
        }
        else if (std::any_of(
                   nu.begin(), nu.end(), [](double x) { return not(x == 0.0 or x > 1.0); }))
          throw std::logic_error("Invalid total fission neutron yield value encountered.\n"
                                 "Only values strictly greater than one, or zero, are "
                                 "permitted.");

        // compute prompt/delayed nu, if needed
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
      } // if nu

      if (fw == "NU_PROMPT_BEGIN")
      {
        Read1DData("NU_PROMPT", nu_prompt, num_groups_, f, ls, ln);

        if (not HasNonZero(nu_prompt))
        {
          Chi::log.Log0Warning() << "The prompt fission neutron yield specified in "
                                 << "\"" << file_name << "\" is uniformly zero..."
                                 << "Clearing it.";
          nu_prompt.clear();
        }
        else if (std::any_of(nu_prompt.begin(),
                             nu_prompt.end(),
                             [](double x) { return not(x == 0.0 or x > 1.0); }))
          throw std::logic_error("Invalid prompt fission neutron yield value encountered.\n"
                                 "Only values strictly greater than one, or zero, are "
                                 "permitted.");
      }

      if (fw == "NU_DELAYED_BEGIN")
      {
        Read1DData("NU_DELAYED", nu_delayed, num_groups_, f, ls, ln);

        if (not HasNonZero(nu_delayed))
        {
          Chi::log.Log0Warning() << "The delayed fission neutron yield specified in "
                                 << "\"" << file_name << "\" is uniformly zero..."
                                 << "Clearing it.";
          nu_prompt.clear();
        }
        else if (not IsNonNegative(nu_delayed))
          throw std::logic_error("Invalid delayed fission neutron yield value encountered.\n"
                                 "Only non-negative values are permitted.");
      }

      if (fw == "BETA_BEGIN")
      {
        Read1DData("BETA", beta, num_groups_, f, ls, ln);

        if (not HasNonZero(beta))
        {
          Chi::log.Log0Warning() << "The delayed neutron fraction specified in "
                                 << "\"" << file_name << "\" is uniformly zero..."
                                 << "Clearing it.";
          beta.clear();
        }
        else if (std::any_of(
                   beta.begin(), beta.end(), [](double x) { return not(x >= 0.0 and x <= 1.0); }))
          throw std::logic_error("Invalid delayed neutron fraction value encountered.\n"
                                 "Only values in the range [0.0, 1.0] are permitted.");

        // compute prompt/delayed nu, if needed
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

      // Fission/Emission Spectra

      if (fw == "CHI_BEGIN")
      {
        Read1DData("CHI", chi, num_groups_, f, ls, ln);

        if (not HasNonZero(chi))
          throw std::logic_error("Invalid steady-state fission spectrum encountered.\n"
                                 "No non-zero values found.");
        if (not IsNonNegative(chi))
          throw std::logic_error("Invalid steady-state fission spectrum value encountered.\n"
                                 "Only non-negative values are permitted.");

        // normalizing
        double sum = std::accumulate(chi.begin(), chi.end(), 0.0);
        std::transform(chi.begin(), chi.end(), chi.begin(), [sum](double& x) { return x / sum; });
      } // if chi

      if (fw == "CHI_PROMPT_BEGIN")
      {
        Read1DData("CHI_PROMPT", chi_prompt, num_groups_, f, ls, ln);

        if (not HasNonZero(chi_prompt))
          throw std::logic_error("Invalid prompt fission spectrum encountered.\n"
                                 "No non-zero values found.");
        if (not IsNonNegative(chi_prompt))
          throw std::logic_error("Invalid prompt fission spectrum value encountered.\n"
                                 "Only non-negative values are permitted.");

        // normalizing
        double sum = std::accumulate(chi_prompt.begin(), chi_prompt.end(), 0.0);
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

        for (unsigned int j = 0; j < num_precursors_; ++j)
        {
          if (not HasNonZero(emission_spectra[j]))
            throw std::logic_error("Invalid delayed emission spectrum encountered for "
                                   "precursor species " +
                                   std::to_string(j) + ".\n" + "No non-zero values found.");
          if (not IsNonNegative(emission_spectra[j]))
            throw std::logic_error("Invalid delayed emission spectrum value encountered "
                                   "for precursor species " +
                                   std::to_string(j) + ".\n" +
                                   "Only non-negative values are permitted.");

          // normalizing
          double sum = std::accumulate(emission_spectra[j].begin(), emission_spectra[j].end(), 0.0);
          std::transform(emission_spectra[j].begin(),
                         emission_spectra[j].end(),
                         emission_spectra[j].begin(),
                         [sum](double& x) { return x / sum; });
        }
      } // if delayed chi

      // Delayed Neutron Precursor Data

      if (num_precursors_ > 0)
      {
        if (fw == "PRECURSOR_DECAY_CONSTANTS_BEGIN")
        {
          Read1DData("PRECURSOR_DECAY_CONSTANTS", decay_constants, num_precursors_, f, ls, ln);

          if (not IsPositive(decay_constants))
            throw std::logic_error("Invalid precursor decay constant value encountered.\n"
                                   "Only strictly positive values are permitted.");
        } // if decay constants

        if (fw == "PRECURSOR_FRACTIONAL_YIELDS_BEGIN")
        {
          Read1DData("PRECURSOR_FRACTIONAL_YIELDS", fractional_yields, num_precursors_, f, ls, ln);

          if (not HasNonZero(fractional_yields))
            throw std::logic_error("Invalid precursor fractional yields encountered.\n"
                                   "No non-zero values found.");
          if (std::any_of(fractional_yields.begin(),
                          fractional_yields.end(),
                          [](double x) { return not(x >= 0.0 and x <= 1.0); }))
            throw std::logic_error("Invalid precursor fractional yield value encountered.\n"
                                   "Only values in the range [0.0, 1.0] are permitted.");

          // normalizing
          double sum = std::accumulate(fractional_yields.begin(), fractional_yields.end(), 0.0);
          std::transform(fractional_yields.begin(),
                         fractional_yields.end(),
                         fractional_yields.begin(),
                         [sum](double& x) { return x / sum; });
        }
      }

      // Transfer Data

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
      throw std::runtime_error("Unknown error encountered in " + std::string(__FUNCTION__));
    }

    word = "";
  } // while not EOF, read each lines
  file.close();

  if (sigma_a_.empty()) ComputeAbsorption();
  ComputeDiffusionParameters();

  // Compute and check fission data

  // determine if the material is fissionable
  is_fissionable_ =
    not sigma_f_.empty() or not nu_sigma_f_.empty() or not production_matrix_.empty();

  // clear fission data if not fissionable
  if (not is_fissionable_)
  {
    sigma_f_.clear();
    nu_sigma_f_.clear();
    nu_prompt_sigma_f_.clear();
    nu_delayed_sigma_f_.clear();
    precursors_.clear();
  } // if not fissionable

  // otherwise, check and set the fission data
  else
  {
    // check vector data inputs
    if (production_matrix_.empty())
    {
      // check for non-delayed fission neutron yield data
      if (nu.empty() and nu_prompt.empty())
        throw std::logic_error("Invalid fission neutron yield specification encountered.\n"
                               "Either the total or prompt fission neutron yield must be "
                               "specified.");

      if (not nu.empty() and not nu_prompt.empty())
        throw std::logic_error("Ambiguous fission neutron yield data specified.\n"
                               "Either the total or prompt fission neutron yield "
                               "must be specified, not both.");

      // check for fission spectrum data
      if (chi.empty() and chi_prompt.empty())
        throw std::logic_error("Invalid fission spectrum specification encountered.\n"
                               "Either the steady-state or prompt fission spectrum must be "
                               "specified.");

      if (not chi.empty() and not chi_prompt.empty())
        throw std::logic_error("Ambiguous fission spectrum data specified.\n"
                               "Either the steady-state or prompt fission spectrum "
                               "must be specified, not both.");

      // check for compatibility
      if ((not nu.empty() and chi.empty()) or (nu.empty() and not chi.empty()) or
          (not nu_prompt.empty() and chi_prompt.empty()) or
          (nu_prompt.empty() and not chi_prompt.empty()))
        throw std::logic_error("Ambiguous fission data specified.\n"
                               "Either the total fission neutron yield and steady-state "
                               "fission spectrum or the prompt fission neutron yield and "
                               "prompt fission spectrum must be specified.");

      // initialize total fission neutron yield
      // do this only when prompt is specified
      if (not nu_prompt.empty())
      {
        nu.assign(num_groups_, 0.0);
        for (unsigned int g = 0; g < num_groups_; ++g)
          nu[g] = nu_prompt[g];
      }

      // check delayed neutron data
      if (num_precursors_ > 0)
      {
        // check that prompt data was specified
        if (chi_prompt.empty() or nu_prompt.empty())
          throw std::logic_error("Invalid fission data specification encountered.\n"
                                 "When delayed neutron precursors are present, "
                                 "prompt fission data must be specified.");

        // check that delayed fission neutron yield was specified
        if (nu_delayed.empty() or decay_constants.empty() or fractional_yields.empty() or
            std::any_of(emission_spectra.begin(),
                        emission_spectra.end(),
                        [](const std::vector<double>& x) { return x.empty(); }))
          throw std::logic_error("Invalid fission data specification encountered.\n"
                                 "When delayed neutron precursors are present, "
                                 "the delayed fission neutron yield, delayed neutron "
                                 "precursor decay constants, fractional yields, and "
                                 "emission spectra must all be specified.");

        // add delayed fission neutron yield to total
        for (unsigned int g = 0; g < num_groups_; ++g)
          nu[g] += nu_delayed[g];

        // add data to precursor structs
        for (unsigned int j = 0; j < num_precursors_; ++j)
        {
          precursors_[j].decay_constant = decay_constants[j];
          precursors_[j].fractional_yield = fractional_yields[j];
          precursors_[j].emission_spectrum = emission_spectra[j];
        }
      }

      // compute fission cross section
      if (sigma_f_.empty() and not nu_sigma_f_.empty())
      {
        sigma_f_.assign(num_groups_, 0.0);
        for (unsigned int g = 0; g < num_groups_; ++g)
          if (nu_sigma_f_[g] > 0.0) sigma_f_[g] = nu_sigma_f_[g] / nu[g];
      }

      // compute total production cross section
      nu_sigma_f_.assign(num_groups_, 0.0);
      for (unsigned int g = 0; g < num_groups_; ++g)
        nu_sigma_f_[g] = nu[g] * sigma_f_[g];

      // compute prompt production cross section
      if (not nu_prompt.empty())
      {
        nu_prompt_sigma_f_.assign(num_groups_, 0.0);
        for (unsigned int g = 0; g < num_groups_; ++g)
          nu_prompt_sigma_f_[g] = nu_prompt[g] * sigma_f_[g];
      }

      // compute delayed production cross section
      if (not nu_delayed.empty())
      {
        nu_delayed_sigma_f_.assign(num_groups_, 0.0);
        for (unsigned int g = 0; g < num_groups_; ++g)
          nu_delayed_sigma_f_[g] = nu_delayed[g] * sigma_f_[g];
      }

      // compute production matrix
      const auto chi_ = not chi_prompt.empty() ? chi_prompt : chi;
      const auto nu_sigma_f = not nu_prompt.empty() ? nu_prompt_sigma_f_ : nu_sigma_f_;

      for (unsigned int g = 0; g < num_groups_; ++g)
      {
        std::vector<double> prod;
        for (unsigned int gp = 0.0; gp < num_groups_; ++gp)
          prod.push_back(chi_[g] * nu_sigma_f[gp]);
        production_matrix_.push_back(prod);
      }
    } // if production_matrix empty

    else
    {
      // TODO: Develop an implementation for multi-particle delayed
      //       neutron data. The primary challenge in this is that
      //       different precursor species exist for neutron-induced
      //       fission than for photo-fission.

      // throw error if precursors are present
      if (num_precursors_ > 0)
        throw std::runtime_error("This setup has not been implemented.\n"
                                 "Currently, when a production matrix is specified, "
                                 "no delayed neutron precursors are allowed.");

      // check for fission cross sections
      if (sigma_f_.empty())
        throw std::logic_error("Invalid fission data specification encountered.\n"
                               "When a production matrix is specified, the fission "
                               "cross sections must also be specified.");

      // compute production cross section
      nu_sigma_f_.assign(num_groups_, 0.0);
      for (unsigned int g = 0; g < num_groups_; ++g)
        for (unsigned int gp = 0; gp < num_groups_; ++gp)
          nu_sigma_f_[gp] += production_matrix_[g][gp];

      // check for reasonable fission neutron yield
      nu.assign(num_groups_, 0.0);
      for (unsigned int g = 0; g < num_groups_; ++g)
        if (sigma_f_[g] > 0.0) nu[g] = nu_sigma_f_[g] / sigma_f_[g];

      if (std::any_of(
            nu.begin(), nu.end(), [](double x) { return not(x == 0.0 or (x > 1.0 and x < 10.0)); }))
        throw std::logic_error("Incompatible fission data encountered.\n"
                               "The estimated fission neutron yield must be either zero "
                               "or in the range (1.0, 10.0).");
    }

    ChiLogicalErrorIf(sigma_f_.empty(),
                      "After reading xs, a fissionable "
                      "material's sigma_f is not defined");
  } // if fissionable
}

void
SingleStateMGXS::ComputeAbsorption()
{
  sigma_a_.assign(num_groups_, 0.0);

  // compute for a pure absorber
  if (transfer_matrices_.empty())
    for (size_t g = 0; g < num_groups_; ++g)
      sigma_a_[g] = sigma_t_[g];

  // estimate from a transfer matrix
  else
  {
    Chi::log.Log0Warning() << "Estimating absorption from the transfer matrices.";

    const auto& S0 = transfer_matrices_[0];
    for (size_t g = 0; g < num_groups_; ++g)
    {
      // estimate the scattering cross section
      double sig_s = 0.0;
      for (size_t row = 0; row < S0.NumRows(); ++row)
      {
        const auto& cols = S0.rowI_indices_[row];
        const auto& vals = S0.rowI_values_[row];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sig_s += vals[t];
            break;
          }
      }

      sigma_a_[g] = sigma_t_[g] - sig_s;

      // TODO: Should negative absorption be allowed?
      if (sigma_a_[g] < 0.0)
        Chi::log.Log0Warning() << "Negative absorption cross section encountered "
                               << "in group " << g << " when estimating from the "
                               << "transfer matrices";
    } // for g
  }   // if scattering present
}

void
SingleStateMGXS::ComputeDiffusionParameters()
{
  if (diffusion_initialized_) return;

  // initialize diffusion data
  diffusion_coeff_.resize(num_groups_, 1.0);
  sigma_s_gtog_.resize(num_groups_, 0.0);
  sigma_removal_.resize(num_groups_, 0.1);

  // perfom computations group-wise
  const auto& S = transfer_matrices_;
  for (unsigned int g = 0; g < num_groups_; ++g)
  {
    // Determine transport correction

    double sig_1 = 0.0;
    if (S.size() > 1)
    {
      for (unsigned int gp = 0; gp < num_groups_; ++gp)
      {
        const auto& cols = S[1].rowI_indices_[gp];
        const auto& vals = S[1].rowI_values_[gp];
        for (size_t t = 0; t < cols.size(); ++t)
          if (cols[t] == g)
          {
            sig_1 += vals[t];
            break;
          }
      } // for gp
    }   // if moment 1 available

    // Compute diffusion coefficient

    if (sig_1 >= sigma_t_[g])
    {
      sig_1 = 0.0;
      Chi::log.Log0Warning() << "Transport corrected diffusion coefficient failed for group " << g
                             << " in call to " << __FUNCTION__ << ". "
                             << "sigma_t=" << sigma_t_[g] << " sigs_g_(m=1)=" << sig_1
                             << ". Setting sigs_g_(m=1) to zero for this group.";
    }

    // compute the diffusion coefficient
    // cap the value for when sig_t - sig_1 is near zero
    diffusion_coeff_[g] = std::fmin(1.0e12, 1.0 / 3.0 / (sigma_t_[g] - sig_1));

    // Determine within group scattering

    if (!S.empty())
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

    sigma_removal_[g] = std::max(0.0, sigma_t_[g] - sigma_s_gtog_[g]);
  } // for g

  diffusion_initialized_ = true;
}

} // namespace chi_physics
