// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/sparse_matrix/sparse_matrix.h"

namespace opensn
{

/// Multi-group cross section container.
class MultiGroupXS : public std::enable_shared_from_this<MultiGroupXS>
{
public:
  /// Default constructor.
  MultiGroupXS()
    : num_groups_(0),
      scattering_order_(0),
      num_precursors_(0),
      is_fissionable_(false),
      adjoint_(false),
      scaling_factor_(1.0),
      diffusion_initialized_(false)
  {
  }

  /**
   * Makes a simple material with a 1-group cross-section set.
   *
   * \param sigma_t Total cross section.
   * \param c Scattering ratio.
   */
  void Initialize(double sigma_t, double c);

  /**
   * Populates the cross section from a combination of others.
   *
   * \param combinations Pairs of cross-section handles and scaling factors.
   */
  void Initialize(std::vector<std::pair<int, double>>& combinations);

  /**
   * Populates transport cross sections from an OpenSn cross-section file.
   *
   * \param file_name Path to the OpenSn cross section file.
   */
  void Initialize(const std::string& file_name);

  /**
   * Populates transport cross sections from an OpenMC cross-section file.
   *
   * \param file_name Path to the HDF5 file.
   * \param dataset_name Dataset inside the file.
   * \param temperature Temperature in Kelvin.
   */
  void
  Initialize(const std::string& file_name, const std::string& dataset_name, double temperature);

  /// Data for a delayed neutron precursor.
  struct Precursor
  {
    double decay_constant = 0.0;
    double fractional_yield = 0.0;
    std::vector<double> emission_spectrum;
  };

  /**
   * Scale the cross sections by the specified factor.
   *
   * @note Scaling factors do not compound. Each time this routine is called, the cross sections
   *       are scaled by the ratio of the argument and the existing scaling factor.
   */
  void SetScalingFactor(const double factor);

  /**
   * Exports the cross-section information to OpenSn format.
   *
   * \param file_name The name of the file to save the cross sections to.
   * \param fission_scaling A factor to scale fission data to. This is generally equal to \f$
   *        1/k_{eff} \f$. Generally, this is done to create exactly critical cross sections for a
   *        transient initial condition.
   */
  void ExportToOpenSnXSFile(const std::string& file_name, const double fission_scaling = 1.0) const;

  /// Returns the number of energy groups.
  size_t GetNumGroups() const { return num_groups_; }

  /// Returns the scattering order.
  size_t GetScatteringOrder() const { return scattering_order_; }

  /// Returns the number of precursor families.
  size_t GetNumPrecursors() const { return num_precursors_; }

  /// Indicates if the material is fissionable.
  bool IsFissionable() const { return is_fissionable_; }

  /// Enable or disable adjoint mode.
  void SetAdjointMode(bool val)
  {
    adjoint_ = val;
    if (adjoint_ and transposed_transfer_matrices_.empty())
      TransposeTransferAndProduction();
  }

  /// Returns true if adjoint mode is enabled.
  bool GetAdjointMode() const { return adjoint_; }

  /// Current scaling factor applied to the data.
  double GetScalingFactor() const { return scaling_factor_; }

  /// Total cross section per group.
  const std::vector<double>& GetSigmaTotal() const { return sigma_t_; }

  /// Absorption cross section per group.
  const std::vector<double>& GetSigmaAbsorption() const { return sigma_a_; }

  /// Transfer matrices ordered by scattering moment.
  const std::vector<SparseMatrix>& GetTransferMatrices() const
  {
    return adjoint_ ? transposed_transfer_matrices_ : transfer_matrices_;
  }

  /// Transfer matrix for a specific scattering moment.
  const SparseMatrix& GetTransferMatrix(unsigned int ell) const
  {
    return adjoint_ ? transposed_transfer_matrices_.at(ell) : transfer_matrices_.at(ell);
  }

  /// Fission spectrum.
  const std::vector<double>& GetChi() const { return chi_; }

  /// Fission cross section per group.
  const std::vector<double>& GetSigmaFission() const { return sigma_f_; }

  /// Neutron production cross section.
  const std::vector<double>& GetNuSigmaF() const { return nu_sigma_f_; }

  /// Prompt neutron production cross section.
  const std::vector<double>& GetNuPromptSigmaF() const { return nu_prompt_sigma_f_; }

  /// Delayed neutron production cross section.
  const std::vector<double>& GetNuDelayedSigmaF() const { return nu_delayed_sigma_f_; }

  /// Total neutron production matrix.
  const std::vector<std::vector<double>>& GetProductionMatrix() const
  {
    return adjoint_ ? transposed_production_matrix_ : production_matrix_;
  }

  /// Delayed neutron precursor data.
  const std::vector<Precursor>& GetPrecursors() const { return precursors_; }

  /// Inverse neutron velocities.
  const std::vector<double>& GetInverseVelocity() const { return inv_velocity_; }

  /// Indicates if diffusion data have been initialized.
  bool DiffusionInitialized() const { return diffusion_initialized_; }

  /// Transport cross section per group.
  const std::vector<double>& GetSigmaTransport() const { return sigma_tr_; }

  /// Diffusion coefficient per group.
  const std::vector<double>& GetDiffusionCoefficient() const { return diffusion_coeff_; }

  /// Removal cross section per group.
  const std::vector<double>& GetSigmaRemoval() const { return sigma_r_; }

  /// Within-group scattering cross section per group.
  const std::vector<double>& GetSigmaSGtoG() const { return sigma_s_gtog_; }

private:
  /// Total number of groups
  size_t num_groups_;
  /// Legendre scattering order
  size_t scattering_order_;
  /// Number of precursors
  size_t num_precursors_;
  /// Is fissionable?
  bool is_fissionable_;
  /// Can be used for adjoint calculations
  bool adjoint_;
  /// An arbitrary scaling factor
  double scaling_factor_ = 1.0;
  /// Evaluation temperature
  double temperature_ = 294.0;
  /// Energy bin boundaries in MeV
  std::vector<double> e_bounds_;
  /// Total cross section
  std::vector<double> sigma_t_;
  /// Absorption cross section
  std::vector<double> sigma_a_;
  /// Fission cross section
  std::vector<double> sigma_f_;
  /// Neutron production due to fission
  std::vector<double> nu_sigma_f_;
  /// Neutron fission spectrum
  std::vector<double> chi_;
  /// Prompt neutron production due to fission
  std::vector<double> nu_prompt_sigma_f_;
  /// Delayed neutron production due to fission
  std::vector<double> nu_delayed_sigma_f_;
  /// Inverse velocity
  std::vector<double> inv_velocity_;
  std::vector<Precursor> precursors_;
  /// Sparse scattering matrix
  std::vector<SparseMatrix> transfer_matrices_;
  std::vector<SparseMatrix> transposed_transfer_matrices_;
  /// Total neutron production matrix
  std::vector<std::vector<double>> production_matrix_;
  std::vector<std::vector<double>> transposed_production_matrix_;

  // Diffusion quantities
  bool diffusion_initialized_;
  /// Transport cross section
  std::vector<double> sigma_tr_;
  /// Transport corrected diffusion coefficient
  std::vector<double> diffusion_coeff_;
  /// Removal cross section
  std::vector<double> sigma_r_;
  /// Within-group scattering cross section
  std::vector<double> sigma_s_gtog_;

  /// Resets all data to the uninitialized state.
  void Reset();

  /// Computes absorption from the total and scattering matrices.
  void ComputeAbsorption();

  /// Computes various diffusion-related quantities.
  void ComputeDiffusionParameters();

  /// Populates transposed transfer and production matrices.
  void TransposeTransferAndProduction();

  /// Check vector for all non-negative values
  bool IsNonNegative(const std::vector<double>& vec)
  {
    return not std::any_of(vec.begin(), vec.end(), [](double x) { return x < 0.0; });
  }

  /// Check vector for all strictly positive values
  bool IsPositive(const std::vector<double>& vec)
  {
    return not std::any_of(vec.begin(), vec.end(), [](double x) { return x <= 0.0; });
  }

  /// Check vector for any non-zero values
  bool HasNonZero(const std::vector<double>& vec)
  {
    return std::any_of(vec.begin(), vec.end(), [](double x) { return x > 0.0; });
  }
};

} // namespace opensn
