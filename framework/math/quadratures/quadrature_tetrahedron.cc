#include "framework/math/quadratures/quadrature_tetrahedron.h"

#include "framework/math/quadratures/conical/quadrature_conical.h"

#include <cmath>
#include <stdexcept>
#include <cassert>

namespace opensn
{

QuadratureTetrahedron::QuadratureTetrahedron(QuadratureOrder order) : Quadrature(order)
{
  double x = 0.0, y = 0.0, z = 0.0;

  switch (order)
  {
    case QuadratureOrder::FIRST:
    {
      x = 0.25;
      y = 0.25;
      z = 0.25;
      qpoints_.emplace_back(x, y, z);
      weights_.push_back(1.0 / 6.0);
      break;
    }
    case QuadratureOrder::SECOND:
    {
      double a = (5.0 + 3.0 * sqrt(5)) / 20.0;
      double b = (5.0 - sqrt(5)) / 20.0;

      x = a;
      y = b;
      z = b;
      qpoints_.emplace_back(x, y, z);
      weights_.push_back(1.0 / 24.0);

      x = b;
      y = a;
      z = b;
      qpoints_.emplace_back(x, y, z);
      weights_.push_back(1.0 / 24.0);

      x = b;
      y = b;
      z = a;
      qpoints_.emplace_back(x, y, z);
      weights_.push_back(1.0 / 24.0);

      x = b;
      y = b;
      z = b;
      qpoints_.emplace_back(x, y, z);
      weights_.push_back(1.0 / 24.0);
      break;
    }
    case QuadratureOrder::THIRD:
    {
      QuadratureConical conical(order);
      conical.Initialize_Conical_Product_Tet();
      qpoints_.swap(conical.qpoints_);
      weights_.swap(conical.weights_);
      break;
    }
    case QuadratureOrder::FOURTH:

    // Walkington's fifth-order 14-point rule from
    // "Quadrature on Simplices of Arbitrary Dimension"
    case QuadratureOrder::FIFTH:
    {
      qpoints_.resize(14);
      weights_.resize(14);

      // permutations of these points and suitably-modified versions of
      // these points are the quadrature point locations
      const double a[3] = {0.31088591926330060980,   // a1 from the paper
                           0.092735250310891226402,  // a2 from the paper
                           0.045503704125649649492}; // a3 from the paper

      // weights.  a[] and wt[] are the only floating-point inputs required
      // for this rule.
      const double wt[3] = {0.018781320953002641800,   // w1 from the paper
                            0.012248840519393658257,   // w2 from the paper
                            0.0070910034628469110730}; // w3 from the paper

      // The first two sets of 4 points are formed in a similar manner
      for (unsigned int i = 0; i < 2; ++i)
      {
        // Where we will insert values into qpoints and weights
        const unsigned int offset = 4 * i;

        // Stuff points and weights values into their arrays
        const double b = 1. - 3. * a[i];

        // Here are the permutations.  Order of these is not important,
        // all have the same weight
        qpoints_[offset + 0] = Vector3(a[i], a[i], a[i]);
        qpoints_[offset + 1] = Vector3(a[i], b, a[i]);
        qpoints_[offset + 2] = Vector3(b, a[i], a[i]);
        qpoints_[offset + 3] = Vector3(a[i], a[i], b);

        // These 4 points all have the same weights
        for (unsigned int j = 0; j < 4; ++j)
          weights_[offset + j] = wt[i];
      } // end for

      {
        // The third set contains 6 points and is formed a little differently
        const unsigned int offset = 8;
        const double b = 0.5 * (1. - 2. * a[2]);

        // Here are the permutations.  Order of these is not important,
        // all have the same weight
        qpoints_[offset + 0] = Vector3(b, b, a[2]);
        qpoints_[offset + 1] = Vector3(b, a[2], a[2]);
        qpoints_[offset + 2] = Vector3(a[2], a[2], b);
        qpoints_[offset + 3] = Vector3(a[2], b, a[2]);
        qpoints_[offset + 4] = Vector3(b, a[2], b);
        qpoints_[offset + 5] = Vector3(a[2], b, b);

        // These 6 points all have the same weights
        for (unsigned int j = 0; j < 6; ++j)
          weights_[offset + j] = wt[2];
      }
      break;
    }

      // This rule is originally from Keast:
      //    Patrick Keast,
      //    Moderate Degree Tetrahedral Quadrature Formulas,
      //    Computer Methods in Applied Mechanics and Engineering,
      //    Volume 55, Number 3, May 1986, pages 339-348.
      //
      // It is accurate on 6th-degree polynomials and has 24 points
      // vs. 64 for the comparable conical product rule.
      //
      // Values copied 24th June 2008 from:
      // http://people.scs.fsu.edu/~burkardt/f_src/keast/keast.f90
    case QuadratureOrder::SIXTH:
    {
      qpoints_.resize(24);
      weights_.resize(24);

      // The raw data for the quadrature rule.
      const std::vector<std::vector<double>> rule_data = {
        {0.356191386222544953e+00, 0.214602871259151684e+00, 0., 0.00665379170969464506e+00},  // 4
        {0.877978124396165982e+00, 0.0406739585346113397e+00, 0., 0.00167953517588677620e+00}, // 4
        {0.0329863295731730594e+00, 0.322337890142275646e+00, 0., 0.00922619692394239843e+00}, // 4
        {0.0636610018750175299e+00,
         0.269672331458315867e+00,
         0.603005664791649076e+00,
         0.00803571428571428248e+00} // 12
      };

      // Now call the keast routine to generate _points and _weights
      KeastRule(rule_data, 4);
      break;
    }
    default:
    {
      QuadratureConical conical(order);
      conical.Initialize_Conical_Product_Tet();
      qpoints_.swap(conical.qpoints_);
      weights_.swap(conical.weights_);
    }
  } // switch order
}

void
QuadratureTetrahedron::KeastRule(const std::vector<std::vector<double>>& rule_data,
                                 const unsigned int n_pts)
{
  auto& _points = qpoints_;
  auto& _weights = weights_;

  typedef Vector3 Point;

  // Like the Dunavant rule, the input data should have 4 columns.  These
  // columns have the following format and implied permutations (w=weight).
  // {a, 0, 0, w} = 1-permutation  (a,a,a)
  // {a, b, 0, w} = 4-permutation  (a,b,b), (b,a,b), (b,b,a), (b,b,b)
  // {a, 0, b, w} = 6-permutation  (a,a,b), (a,b,b), (b,b,a), (b,a,b), (b,a,a), (a,b,a)
  // {a, b, c, w} = 12-permutation (a,a,b), (a,a,c), (b,a,a), (c,a,a), (a,b,a), (a,c,a)
  //                               (a,b,c), (a,c,b), (b,a,c), (b,c,a), (c,a,b), (c,b,a)

  // Always insert into the points & weights vector relative to the offset
  unsigned int offset = 0;

  for (unsigned int p = 0; p < n_pts; ++p)
  {

    // There must always be a non-zero entry to start the row
    assert(rule_data[p][0] != static_cast<double>(0.0));

    // A zero weight may imply you did not set up the raw data correctly
    assert(rule_data[p][3] != static_cast<double>(0.0));

    // What kind of point is this?
    // One non-zero entry in first 3 cols   ? 1-perm (centroid) point = 1
    // Two non-zero entries in first 3 cols ? 3-perm point            = 3
    // Three non-zero entries               ? 6-perm point            = 6
    unsigned int pointtype = 1;

    if (rule_data[p][1] != static_cast<double>(0.0))
    {
      if (rule_data[p][2] != static_cast<double>(0.0))
        pointtype = 12;
      else
        pointtype = 4;
    }
    else
    {
      // The second entry is zero.  What about the third?
      if (rule_data[p][2] != static_cast<double>(0.0))
        pointtype = 6;
    }

    switch (pointtype)
    {
      case 1:
      {
        // Be sure we have enough space to insert this point
        assert(offset + 0 < _points.size());

        const double a = rule_data[p][0];

        // The point has only a single permutation (the centroid!)
        _points[offset + 0] = Point(a, a, a);

        // The weight is always the last entry in the row.
        _weights[offset + 0] = rule_data[p][3];

        offset += pointtype;
        break;
      }

      case 4:
      {
        // Be sure we have enough space to insert these points
        assert(offset + 3 < _points.size());

        const double a = rule_data[p][0];
        const double b = rule_data[p][1];
        const double wt = rule_data[p][3];

        // Here it's understood the second entry is to be used twice, and
        // thus there are three possible permutations.
        _points[offset + 0] = Point(a, b, b);
        _points[offset + 1] = Point(b, a, b);
        _points[offset + 2] = Point(b, b, a);
        _points[offset + 3] = Point(b, b, b);

        for (unsigned int j = 0; j < pointtype; ++j)
          _weights[offset + j] = wt;

        offset += pointtype;
        break;
      }

      case 6:
      {
        // Be sure we have enough space to insert these points
        assert(offset + 5 < _points.size());

        const double a = rule_data[p][0];
        const double b = rule_data[p][2];
        const double wt = rule_data[p][3];

        // Three individual entries with six permutations.
        _points[offset + 0] = Point(a, a, b);
        _points[offset + 1] = Point(a, b, b);
        _points[offset + 2] = Point(b, b, a);
        _points[offset + 3] = Point(b, a, b);
        _points[offset + 4] = Point(b, a, a);
        _points[offset + 5] = Point(a, b, a);

        for (unsigned int j = 0; j < pointtype; ++j)
          _weights[offset + j] = wt;

        offset += pointtype;
        break;
      }

      case 12:
      {
        // Be sure we have enough space to insert these points
        assert(offset + 11 < _points.size());

        const double a = rule_data[p][0];
        const double b = rule_data[p][1];
        const double c = rule_data[p][2];
        const double wt = rule_data[p][3];

        // Three individual entries with six permutations.
        _points[offset + 0] = Point(a, a, b);
        _points[offset + 6] = Point(a, b, c);
        _points[offset + 1] = Point(a, a, c);
        _points[offset + 7] = Point(a, c, b);
        _points[offset + 2] = Point(b, a, a);
        _points[offset + 8] = Point(b, a, c);
        _points[offset + 3] = Point(c, a, a);
        _points[offset + 9] = Point(b, c, a);
        _points[offset + 4] = Point(a, b, a);
        _points[offset + 10] = Point(c, a, b);
        _points[offset + 5] = Point(a, c, a);
        _points[offset + 11] = Point(c, b, a);

        for (unsigned int j = 0; j < pointtype; ++j)
          _weights[offset + j] = wt;

        offset += pointtype;
        break;
      }

      default:
        throw std::invalid_argument(std::string(__FUNCTION__) +
                                    ": Don't know what to do with this many permutation points!");
    }
  }
}

} // namespace opensn
