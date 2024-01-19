#pragma once

#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"

namespace opensn
{

/**
 * A class for handling multi-group cross sections.
 *
 * See doc_ChiFormatXS.h for cross-section format
 */
class SingleStateMGXS : public MultiGroupXS
{
protected:
  typedef std::vector<std::pair<double, double>> AnglePairs;

public:
  SingleStateMGXS()
    : MultiGroupXS(),
      num_groups_(0),
      scattering_order_(0),
      num_precursors_(0),
      diffusion_initialized_(false),
      scattering_initialized_(false)
  {
  }

  /**
   * Makes a simple material with no transfer matrix just sigma_t.
   */
  void MakeSimple0(unsigned int num_groups, double sigma_t);

  /**
   * Makes a simple material with isotropic transfer matrix (L=0) and mostly down scattering but
   * with a few of the last groups subject to up-scattering.
   */
  void MakeSimple1(unsigned int num_groups, double sigma_t, double c);

  /**
   * Populates the cross section from a combination of others.
   */
  void MakeCombined(std::vector<std::pair<int, double>>& combinations);

private:
  void Clear();

public:
  /**
   * This method populates a transport cross section from a Chi cross section file.
   */
  void MakeFromChiXSFile(const std::string& file_name);

private:
  void ComputeAbsorption();
  void ComputeDiffusionParameters();

public:
  // Accessors
  size_t NumGroups() const override { return num_groups_; }
  size_t ScatteringOrder() const override { return scattering_order_; }
  size_t NumPrecursors() const override { return num_precursors_; }
  bool IsFissionable() const override { return is_fissionable_; }

  const std::vector<double>& SigmaTotal() const override { return sigma_t_; }
  const std::vector<double>& SigmaAbsorption() const override { return sigma_a_; }
  const std::vector<SparseMatrix>& TransferMatrices() const override { return transfer_matrices_; }
  const SparseMatrix& TransferMatrix(unsigned int ell) const override
  {
    return transfer_matrices_.at(ell);
  }

  const std::vector<double>& SigmaFission() const override { return sigma_f_; }
  const std::vector<double>& NuSigmaF() const override { return nu_sigma_f_; }
  const std::vector<double>& NuPromptSigmaF() const override { return nu_prompt_sigma_f_; }
  const std::vector<double>& NuDelayedSigmaF() const override { return nu_delayed_sigma_f_; }
  const std::vector<std::vector<double>>& ProductionMatrix() const override
  {
    return production_matrix_;
  }
  const std::vector<Precursor>& Precursors() const override { return precursors_; }

  const std::vector<double>& InverseVelocity() const override { return inv_velocity_; }

  bool DiffusionInitialized() const override { return diffusion_initialized_; }
  const std::vector<double>& SigmaTransport() const override { return sigma_tr_; }
  const std::vector<double>& DiffusionCoefficient() const override { return diffusion_coeff_; }
  const std::vector<double>& SigmaRemoval() const override { return sigma_r_; }
  const std::vector<double>& SigmaSGtoG() const override { return sigma_s_gtog_; }

protected:
  size_t num_groups_ = 0;       ///< Total number of groups
  size_t scattering_order_ = 0; ///< Legendre scattering order
  size_t num_precursors_ = 0;   ///< Number of precursors
  bool is_fissionable_ = false;
  std::vector<std::vector<double>> e_bounds_; ///< Energy bin boundaries in MeV

  std::vector<double> sigma_t_; ///< Total cross section
  std::vector<double> sigma_a_; ///< Absorption cross section
  std::vector<SparseMatrix> transfer_matrices_;

  std::vector<double> sigma_f_; ///< Fission cross section
  std::vector<double> nu_sigma_f_;
  std::vector<double> nu_prompt_sigma_f_;
  std::vector<double> nu_delayed_sigma_f_;
  std::vector<std::vector<double>> production_matrix_;
  std::vector<Precursor> precursors_;

  std::vector<double> inv_velocity_;

  // Diffusion quantities
  bool diffusion_initialized_ = false;
  std::vector<double> sigma_tr_;        ///< Transport cross section
  std::vector<double> diffusion_coeff_; ///< Transport corrected diffusion coefficient
  std::vector<double> sigma_r_;         ///< Removal cross section
  std::vector<double> sigma_s_gtog_;    ///< Within-group scattering cross section

  // Monte-Carlo quantities
private:
  bool scattering_initialized_ = false;
  std::vector<std::vector<double>> cdf_gprime_g_;
  std::vector<std::vector<AnglePairs>> scat_angles_gprime_g_;
};

} // namespace opensn
