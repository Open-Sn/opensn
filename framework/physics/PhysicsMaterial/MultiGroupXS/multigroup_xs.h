#ifndef MULTIGROUP_XS_H
#define MULTIGROUP_XS_H

#include "physics/PhysicsMaterial/material_property_base.h"
#include "math/SparseMatrix/chi_math_sparse_matrix.h"


namespace chi_physics
{

//######################################################################
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

  MultiGroupXS()
      : MaterialProperty(PropertyType::TRANSPORT_XSECTIONS)
  {}

  void ExportToChiXSFile(const std::string& file_name,
                         const double fission_scaling = 1.0) const;
  void PushLuaTable(lua_State* L) const override;

  virtual const unsigned int NumGroups() const = 0;

  virtual const unsigned int ScatteringOrder() const = 0;

  virtual const unsigned int NumPrecursors() const = 0;

  virtual const bool IsFissionable() const = 0;

  virtual const bool DiffusionInitialized() const = 0;

  virtual const bool ScatteringInitialized() const = 0;

  virtual const std::vector<double>& SigmaTotal() const = 0;

  virtual const std::vector<double>& SigmaAbsorption() const = 0;

  virtual const std::vector<double>& SigmaFission() const = 0;

  virtual const std::vector<double>& NuSigmaF() const = 0;

  virtual const std::vector<double>& NuPromptSigmaF() const = 0;

  virtual const std::vector<double>& NuDelayedSigmaF() const = 0;

  virtual const std::vector<double>& InverseVelocity() const = 0;

  virtual const std::vector <chi_math::SparseMatrix>&
  TransferMatrices() const = 0;

  virtual const chi_math::SparseMatrix&
  TransferMatrix(unsigned int ell) const = 0;

  virtual const std::vector <std::vector<double>> ProductionMatrix() const = 0;

  virtual const std::vector <Precursor>& Precursors() const = 0;

  virtual const std::vector<double>& DiffusionCoefficient() const = 0;

  virtual std::vector<double> SigmaTransport() const = 0;

  virtual const std::vector<double>& SigmaRemoval() const = 0;

  virtual const std::vector<double>& SigmaSGtoG() const = 0;
};

}//namespace chi_physics

#endif //MULTIGROUP_XS_H
