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

    // Quadrature-specific order parameter
    // For Product: N in Sn notation
    // For Lebedev: Lebedev order (3, 5, 7, 9, etc.)
    // For SLDFESQ: Uniform Refinement level
    unsigned int quadrature_order = 0;

    // Old Information
    unsigned int N = 0;
    unsigned int num_polar = 0;
    unsigned int num_azimuthal = 0;
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

  // Lebedev quadrature rules
  static std::vector<AngularQuadrature::HarmonicIndices>
  SelectLebedev(const SelectionParameters& params);

  // SLDFESQ quadrature rules
  static std::vector<AngularQuadrature::HarmonicIndices>
  SelectSLDFESQ(const SelectionParameters& params);

  // Helper functions for determining quadrature order from number of angles
  static int DetermineLebedevOrder(size_t num_angles);
  static int DetermineSLDFELevel(size_t num_angles);

  // Standard fallback for non-Galerkin methods
  static std::vector<AngularQuadrature::HarmonicIndices>
  SelectStandard(const SelectionParameters& params);
};

} // namespace opensn