#pragma once

#include "framework/object.h"
#include "framework/mesh/mesh.h"

#include <vector>

namespace chi_math
{
// clang-format off
  enum class QuadratureOrder : int {
    CONSTANT = 0, FIRST = 1, SECOND = 2, THIRD = 3,
    FOURTH = 4, FIFTH = 5, SIXTH = 6, SEVENTH = 7,
    EIGHTH = 8, NINTH = 9, TENTH = 10, ELEVENTH = 11,
    TWELFTH = 12, THIRTEENTH = 13, FOURTEENTH = 14, FIFTEENTH = 15,
    SIXTEENTH = 16, SEVENTEENTH = 17, EIGHTTEENTH = 18, NINETEENTH = 19,
    TWENTIETH = 20, TWENTYFIRST = 21, TWENTYSECOND = 22, TWENTYTHIRD = 23,
    TWENTYFOURTH = 24, TWENTYFIFTH = 25, TWENTYSIXTH = 26, TWENTYSEVENTH = 27,
    TWENTYEIGHTH = 28, TWENTYNINTH = 29, THIRTIETH = 30, THIRTYFIRST = 31,
    THIRTYSECOND = 32, THIRTYTHIRD = 33, THIRTYFOURTH = 34, THIRTYFIFTH = 35,
    THIRTYSIXTH = 36, THIRTYSEVENTH = 37, THIRTYEIGHTH = 38, THIRTYNINTH = 39,
    FORTIETH = 40, FORTYFIRST = 41, FORTYSECOND = 42, FORTYTHIRD = 43,
    INVALID_ORDER
  };
// clang-format on
typedef chi_mesh::Vector3 QuadraturePointXYZ;
class Quadrature;
} // namespace chi_math

/**Parent class for quadratures.*/
class chi_math::Quadrature : public ChiObject
{
public:
  QuadratureOrder order_;
  std::vector<chi_math::QuadraturePointXYZ> qpoints_;
  std::vector<double> weights_;

  static chi::InputParameters GetInputParameters();

protected:
  /**Interval on which the quadrature is defined
   * (relevant for one-dimensional quadratures only).*/
  std::pair<double, double> range_;
  bool verbose_ = false;

protected:
  explicit Quadrature(const chi::InputParameters& params);
  explicit Quadrature(QuadratureOrder in_order) : order_(in_order), range_({0, 0}) {}

public:
  /**Get the range on which the quadrature is defined
   * (relevant for one-dimensional quadratures only).*/
  const std::pair<double, double>& GetRange() const { return range_; }
  /**Set the range on which the quadrature is defined.
   * (relevant for one-dimensional quadratures only).
   * Note that calling this method results in translation
   * of the abscissae and scaling of the weights.*/
  void SetRange(const std::pair<double, double>& in_range);
};
