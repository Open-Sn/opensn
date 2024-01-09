#pragma once

#include "framework/physics/physics_material/material_property_base.h"
#include "framework/math/sparse_matrix/math_sparse_matrix.h"

namespace opensn
{

class MultiGroupXS : public MaterialProperty
{
public:
  /**
   * A struct containing data for a delayed neutron precursor.
   */
  struct Precursor
  {
    double decay_constant = 0.0;
    double fractional_yield = 0.0;
    std::vector<double> emission_spectrum;
  };

  MultiGroupXS() : MaterialProperty(PropertyType::TRANSPORT_XSECTIONS) {}

  /**
   * Exports the cross section information to ChiTech format.
   *
   * \param file_name The name of the file to save the cross sections to.
   * \param fission_scaling A factor to scale fission data to. This is
   *      generally equal to \f$ 1/k_{eff} \f$. Generally, this is done to
   *      create exactly critical cross sections for a transient initial
   *      condition.
   */
  void ExportToChiXSFile(const std::string& file_name, const double fission_scaling = 1.0) const;

  virtual size_t NumGroups() const = 0;
  virtual size_t ScatteringOrder() const = 0;
  virtual size_t NumPrecursors() const = 0;
  virtual bool IsFissionable() const = 0;
  virtual double ScalingFactor() const = 0;

  virtual const std::vector<double>& SigmaTotal() const = 0;
  virtual const std::vector<double>& SigmaAbsorption() const = 0;
  virtual const std::vector<SparseMatrix>& TransferMatrices() const = 0;
  virtual const SparseMatrix& TransferMatrix(unsigned int ell) const = 0;

  virtual const std::vector<double>& SigmaFission() const = 0;
  virtual const std::vector<double>& NuSigmaF() const = 0;
  virtual const std::vector<double>& NuPromptSigmaF() const = 0;
  virtual const std::vector<double>& NuDelayedSigmaF() const = 0;
  virtual const std::vector<std::vector<double>>& ProductionMatrix() const = 0;
  virtual const std::vector<Precursor>& Precursors() const = 0;

  virtual const std::vector<double>& InverseVelocity() const = 0;

  virtual bool DiffusionInitialized() const = 0;
  virtual const std::vector<double>& SigmaTransport() const = 0;
  virtual const std::vector<double>& DiffusionCoefficient() const = 0;
  virtual const std::vector<double>& SigmaRemoval() const = 0;
  virtual const std::vector<double>& SigmaSGtoG() const = 0;
};

} // namespace opensn
