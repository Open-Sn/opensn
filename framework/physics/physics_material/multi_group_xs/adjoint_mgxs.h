#pragma once

#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"

namespace opensn
{

/**
 * A wrapper class for obtaining multi-group cross section data for
 * adjoint simulations.
 *
 * In adjoint simulations, the transfer and production matrices are transposed.
 * While the respective matrices could be queried by simply flipping the
 * indices, access attempts in this fashion are quite costly. Rather, this
 * class precomputes and stores these transpose operators. Along with this, a
 * reference to an instance of MultiGroupXS is stored. In this class, accessors
 * for vector data call the respective accessor from MultiGroupXS and accessors
 * for transfer and production matrices access the respective transposed data
 * stored in this class.
 */
class AdjointMGXS : public MultiGroupXS
{
private:
  const MultiGroupXS& xs_;
  std::vector<SparseMatrix> transposed_transfer_matrices_;
  std::vector<std::vector<double>> transposed_production_matrices_;

public:
  AdjointMGXS() = delete;
  AdjointMGXS(const AdjointMGXS&) = delete;
  AdjointMGXS(AdjointMGXS&&) = delete;

  explicit AdjointMGXS(const MultiGroupXS& xs);

  // Accessors
  const unsigned int NumGroups() const override { return xs_.NumGroups(); }

  const unsigned int ScatteringOrder() const override { return xs_.ScatteringOrder(); }

  const unsigned int NumPrecursors() const override { return xs_.NumPrecursors(); }

  const bool IsFissionable() const override { return xs_.IsFissionable(); }

  const bool DiffusionInitialized() const override { return xs_.DiffusionInitialized(); }

  const bool ScatteringInitialized() const override { return xs_.ScatteringInitialized(); }

  const std::vector<double>& SigmaTotal() const override { return xs_.SigmaTotal(); }

  const std::vector<double>& SigmaAbsorption() const override { return xs_.SigmaAbsorption(); }

  const std::vector<double>& SigmaFission() const override { return xs_.SigmaFission(); }

  const std::vector<double>& NuSigmaF() const override { return xs_.NuSigmaF(); }

  const std::vector<double>& NuPromptSigmaF() const override { return xs_.NuPromptSigmaF(); }

  const std::vector<double>& NuDelayedSigmaF() const override { return xs_.NuDelayedSigmaF(); }

  const std::vector<double>& InverseVelocity() const override { return xs_.InverseVelocity(); }

  const std::vector<SparseMatrix>& TransferMatrices() const override
  {
    return transposed_transfer_matrices_;
  }

  const SparseMatrix& TransferMatrix(unsigned int ell) const override
  {
    return transposed_transfer_matrices_.at(ell);
  }

  const std::vector<std::vector<double>> ProductionMatrix() const override
  {
    return transposed_production_matrices_;
  }

  const std::vector<Precursor>& Precursors() const override { return xs_.Precursors(); }

  const std::vector<double>& DiffusionCoefficient() const override
  {
    return xs_.DiffusionCoefficient();
  }

  std::vector<double> SigmaTransport() const override { return xs_.SigmaTransport(); }

  const std::vector<double>& SigmaRemoval() const override { return xs_.SigmaRemoval(); }

  const std::vector<double>& SigmaSGtoG() const override { return xs_.SigmaSGtoG(); }
};

} // namespace opensn
