// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include <vector>
#include <functional>

namespace opensn
{

/**
 * \brief Rule-based harmonic selection for different quadrature types.
 *
 * This class implements the specific mathematical rules for selecting
 * spherical harmonics based on quadrature type, dimension, and parameters.
 */
class HarmonicSelectionRules
{
public:
  struct SelectionParameters
  {
    AngularQuadratureType quadrature_type;
    OperatorConstructionMethod construction_method;
    unsigned int dimension;
    unsigned int scattering_order;
    size_t num_angles;

    // Cartesian Product quadrature parameters
    unsigned int N = 0;             // "N" in SN Quadratures
    unsigned int num_polar = 0;     // Number of polar angles (optional)
    unsigned int num_azimuthal = 0; // Number of azimuthal angles (optional)
  };

  /**
   * \brief Select harmonics using predetermined mathematical rules
   *
   * \param params Selection parameters including quadrature type and dimensions
   * \return Vector of selected harmonic indices
   */
  static std::vector<AngularQuadrature::HarmonicIndices>
  SelectHarmonics(const SelectionParameters& params);

private:
  // Cartesian Product rule implementations
  static std::vector<AngularQuadrature::HarmonicIndices>
  Select2DCartesianProduct(const SelectionParameters& params);

  static std::vector<AngularQuadrature::HarmonicIndices>
  Select3DCartesianProduct(const SelectionParameters& params);

  // Rule evaluation helpers for Cartesian Products
  static bool CheckCartesian2DRules(unsigned int ell, int m, unsigned int N);

  static bool CheckCartesian3DRules(unsigned int ell, int m, unsigned int N);

  // Standard fallback for non-Galerkin methods
  static std::vector<AngularQuadrature::HarmonicIndices>
  SelectStandard(const SelectionParameters& params);
};

} // namespace opensn