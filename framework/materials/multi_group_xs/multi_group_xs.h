// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/material_property.h"
#include "framework/math/sparse_matrix/sparse_matrix.h"

namespace opensn
{

class MultiGroupXS : public MaterialProperty
{
public:
  MultiGroupXS()
    : MaterialProperty(PropertyType::TRANSPORT_XSECTIONS),
      num_groups_(0),
      scattering_order_(0),
      num_precursors_(0),
      is_fissionable_(false),
      adjoint_(false),
      scaling_factor_(1.0),
      diffusion_initialized_(false)
  {
  }

  /// Makes a simple material with a 1-group cross-section set.
  void Initialize(double sigma_t, double c);

  /// Populates the cross section from a combination of others.
  void Initialize(std::vector<std::pair<int, double>>& combinations);

  /// This method populates transport cross sections from an OpenSn cross-section file.
  void Initialize(const std::string& file_name);

  /// This method populates transport cross sections from an OpenMC cross-section file.
  void
  Initialize(const std::string& file_name, const std::string& dataset_name, double temperature);

  /// A struct containing data for a delayed neutron precursor.
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

  size_t NumGroups() const { return num_groups_; }

  size_t ScatteringOrder() const { return scattering_order_; }

  size_t NumPrecursors() const { return num_precursors_; }

  bool IsFissionable() const { return is_fissionable_; }

  void SetAdjointMode(bool val)
  {
    adjoint_ = val;
    if (adjoint_ and transposed_transfer_matrices_.empty())
      TransposeTransferAndProduction();
  }

  bool GetAdjointMode() const { return adjoint_; }

  double ScalingFactor() const { return scaling_factor_; }

  const std::vector<double>& SigmaTotal() const { return sigma_t_; }

  const std::vector<double>& SigmaAbsorption() const { return sigma_a_; }

  const std::vector<SparseMatrix>& TransferMatrices() const
  {
    return adjoint_ ? transposed_transfer_matrices_ : transfer_matrices_;
  }

  const SparseMatrix& TransferMatrix(unsigned int ell) const
  {
    return adjoint_ ? transposed_transfer_matrices_.at(ell) : transfer_matrices_.at(ell);
  }

  const std::vector<double>& Chi() const { return chi_; }

  const std::vector<double>& SigmaFission() const { return sigma_f_; }

  const std::vector<double>& NuSigmaF() const { return nu_sigma_f_; }

  const std::vector<double>& NuPromptSigmaF() const { return nu_prompt_sigma_f_; }

  const std::vector<double>& NuDelayedSigmaF() const { return nu_delayed_sigma_f_; }

  const std::vector<std::vector<double>>& ProductionMatrix() const
  {
    return adjoint_ ? transposed_production_matrix_ : production_matrix_;
  }

  const std::vector<Precursor>& Precursors() const { return precursors_; }

  const std::vector<double>& InverseVelocity() const { return inv_velocity_; }

  bool DiffusionInitialized() const { return diffusion_initialized_; }

  const std::vector<double>& SigmaTransport() const { return sigma_tr_; }

  const std::vector<double>& DiffusionCoefficient() const { return diffusion_coeff_; }

  const std::vector<double>& SigmaRemoval() const { return sigma_r_; }

  const std::vector<double>& SigmaSGtoG() const { return sigma_s_gtog_; }

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

  void Reset();

  void ComputeAbsorption();

  void ComputeDiffusionParameters();

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
