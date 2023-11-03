#pragma once

#include <map>
#include <vector>

#include "angular_quadrature_base.h"

namespace chi_math
{
enum class ProductQuadratureType
{
  UNKNOWN = 0,
  GAUSS_LEGENDRE = 1,
  GAUSS_CHEBYSHEV = 2,
  GAUSS_LEGENDRE_LEGENDRE = 3,
  GAUSS_LEGENDRE_CHEBYSHEV = 4,
  CUSTOM_QUADRATURE = 5,
};

class ProductQuadrature;
class AngularQuadratureProdGL;
class AngularQuadratureProdGLL;
class AngularQuadratureProdGLC;
class AngularQuadratureProdCustom;
} // namespace chi_math

//######################################################### Class def
/** Class for product quadratures*/
class chi_math::ProductQuadrature : public chi_math::AngularQuadrature
{
public:
  std::vector<double> polar_ang_;
  std::vector<double> azimu_ang_;

protected:
  /** Linear indices of ordered directions mapped to polar level. */
  std::map<unsigned int, std::vector<unsigned int>> map_directions_;

protected:
  ProductQuadrature() : AngularQuadrature(chi_math::AngularQuadratureType::ProductQuadrature) {}

public:
  ~ProductQuadrature() override = default;

  void AssembleCosines(const std::vector<double>& azimuthal,
                       const std::vector<double>& polar,
                       const std::vector<double>& in_weights,
                       bool verbose);

  void OptimizeForPolarSymmetry(double normalization) override;
  /**Obtains the abscissae index given the indices of the
   * polar angle index and the azimuthal angle index.*/
  unsigned int GetAngleNum(const unsigned int polar_angle_index,
                           const unsigned int azimu_angle_index) const
  {
    return map_directions_.at(polar_angle_index)[azimu_angle_index];
  }
  /** Return constant reference to map_directions. */
  const std::map<unsigned int, std::vector<unsigned int>>& GetDirectionMap() const
  {
    return map_directions_;
  }
};

//######################################################### Class def
class chi_math::AngularQuadratureProdGL : public chi_math::ProductQuadrature
{
public:
  explicit AngularQuadratureProdGL(int Np, bool verbose = false);
};

//######################################################### Class def
class chi_math::AngularQuadratureProdGLL : public chi_math::ProductQuadrature
{
public:
  explicit AngularQuadratureProdGLL(int Na, int Np, bool verbose = false);
};

//######################################################### Class def
class chi_math::AngularQuadratureProdGLC : public chi_math::ProductQuadrature
{
public:
  explicit AngularQuadratureProdGLC(int Na, int Np, bool verbose = false);
};

//######################################################### Class def
class chi_math::AngularQuadratureProdCustom : public chi_math::ProductQuadrature
{
public:
  AngularQuadratureProdCustom(const std::vector<double>& azimuthal,
                              const std::vector<double>& polar,
                              const std::vector<double>& in_weights,
                              bool verbose);
};


