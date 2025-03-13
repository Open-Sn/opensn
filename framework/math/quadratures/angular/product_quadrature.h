// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/parameters/input_parameters.h"
#include <map>
#include <vector>

namespace opensn
{

/// Class for product quadratures
class ProductQuadrature : public AngularQuadrature
{
protected:
  /// Linear indices of ordered directions mapped to polar level.
  std::map<unsigned int, std::vector<unsigned int>> map_directions_;

  ProductQuadrature(int dimension)
    : AngularQuadrature(AngularQuadratureType::ProductQuadrature, dimension)
  {
  }

  virtual ~ProductQuadrature() = default;

  /**
   * Optimizes the angular quadrature for polar symmetry by removing all the direction with downward
   * pointing polar angles.
   *
   * \param normalization float. (Optional) The default is a negative number which does not apply
   *        any normalization. If a positive number is provided, the weights will be normalized to
   *        sum to this number.
   */
  void OptimizeForPolarSymmetry(double normalization);

  /// Initializes the quadrature with custom angles and weights.
  void AssembleCosines(const std::vector<double>& azimuthal,
                       const std::vector<double>& polar,
                       const std::vector<double>& weights,
                       bool verbose);

public:
  std::vector<double> polar_ang;
  std::vector<double> azimu_ang;

  /**
   * Obtains the abscissae index given the indices of the polar angle index and the azimuthal angle
   * index.
   */
  unsigned int GetAngleNum(const unsigned int polar_angle_index,
                           const unsigned int azimu_angle_index) const
  {
    return map_directions_.at(polar_angle_index)[azimu_angle_index];
  }

  /// Return constant reference to map_directions.
  const std::map<unsigned int, std::vector<unsigned int>>& GetDirectionMap() const
  {
    return map_directions_;
  }
};

class GLProductQuadrature1DSlab : public ProductQuadrature
{
public:
  /// Constructors for 1D slab Gauss-Legendre product quadrature
  explicit GLProductQuadrature1DSlab(const InputParameters& params);
  GLProductQuadrature1DSlab(int Npolar, bool verbose = false);

  virtual ~GLProductQuadrature1DSlab() = default;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<GLProductQuadrature1DSlab> Create(const ParameterBlock& params);

private:
  void Initialize(int Npolar, bool verbose);
};

class GLCProductQuadrature2DXY : public ProductQuadrature
{
public:
  /// Constructor for 2D XY Gauss-Legendre Chebyshev product quadrature
  explicit GLCProductQuadrature2DXY(const InputParameters& params);
  GLCProductQuadrature2DXY(int Npolar, int Nazimuthal, bool verbose = false);

  virtual ~GLCProductQuadrature2DXY() = default;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<GLCProductQuadrature2DXY> Create(const ParameterBlock& params);

private:
  void Initialize(int Npolar, int Nazimuthal, bool verbose);
};

class GLCProductQuadrature3DXYZ : public ProductQuadrature
{
public:
  /// Constructor for 3D XYZ Gauss-Legendre Chebyshev product quadrature
  explicit GLCProductQuadrature3DXYZ(const InputParameters& params);
  GLCProductQuadrature3DXYZ(int Npolar, int Nazimuthal, bool verbose = false);

  virtual ~GLCProductQuadrature3DXYZ() = default;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<GLCProductQuadrature3DXYZ> Create(const ParameterBlock& params);

private:
  void Initialize(int Npolar, int Nazimuthal, bool verbose);
};

} // namespace opensn
