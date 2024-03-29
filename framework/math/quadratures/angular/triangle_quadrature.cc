#include "framework/math/quadratures/angular/triangle_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/math.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

namespace opensn
{

TriangleQuadrature::TriangleQuadrature(unsigned int method, unsigned int sn, unsigned int moments)
  : method_(method), sn_(sn), moments_(moments)
{
  TriangleInit();
}

void
TriangleQuadrature::TriangleInit()
{
  log.Log0() << "Method: " << method_ << "\nsn " << sn_;

  if (method_ != 1 and method_ != 2 and method_ != 3 and method_ != 0)
  {
    log.Log0Error() << "TriangleQuarature: " << method_ << " is not a valid method.\n"
                    << "Method must be 0, 1, 2, or 3." << std::endl;
    Exit(510);
  }
  if (moments_ < 0 or moments_ > sn_)
  {
    log.Log0Error() << "TriangleQuarature: " << moments_ << " must be > 0 and < " << sn_
                    << std::endl;
    Exit(510);
  }

  // Get the gauss points for the z axis. The number of points for GL is twice
  // that of the quadrature.
  Vector3 new_omega;
  const auto old_omega = GaussLegendreQuadrature(sn_);

  // Formulate the triangular quadrature
  VecDbl newZi, newWeights;
  newZi.reserve(old_omega.qpoints.size());
  newWeights.reserve(old_omega.qpoints.size());

  for (size_t pos = 0; pos < old_omega.qpoints.size(); ++pos)
  {
    if (old_omega.qpoints[pos].x < 0)
      continue;
    newZi.push_back(old_omega.qpoints[pos].x);
    newWeights.push_back(old_omega.weights[pos]);
  }

  int num_div = 1;
  int weightPos = 0;
  for (auto u : old_omega.qpoints)
  {
    double deltaVPhi = M_PI / (2.0 * static_cast<double>(num_div));
    // When the GaussLegendre quadrature gives us the x points we will
    // use for our z positions on the unit sphere, they are ordered in
    // descending order - largest magnitude towards 0 on the negative side,
    // and ascending order on the positive side of 0
    // The weights_ are defined using the weights_ given by the GL quadrature
    // which is the position in the weights_ that weightPos keeps track of.
    // This will ignore the positive x values, and just use the descending
    // order given by the quadrature and use the absolute value of the x values.

    if (u.x >= 0.0)
      break;

    for (int v = 0; v < num_div; ++v)
    {
      double new_z_value = abs(u.x);
      double phi = deltaVPhi / 2.0 + (double)v * deltaVPhi;
      double theta = acos(new_z_value);
      double sinTheta = sqrt(1 - new_z_value * new_z_value);
      double weightCurrent = old_omega.weights[weightPos] / (num_div);

      new_omega.x = sinTheta * cos(phi);
      new_omega.y = sinTheta * sin(phi);
      new_omega.z = new_z_value;
      weights_.push_back(old_omega.weights[weightPos] * (M_PI / (2.0 * num_div)));
      omegas_.emplace_back(new_omega);
      abscissae_.emplace_back(phi, theta);
    }

    weightPos++;
    num_div++;
  }

  // This is the number of points in 1 octant
  size_t octSize = weights_.size();
  for (int octant = 1; octant <= 3; ++octant)
  {
    // This is how much should be added to each phi
    // to get the orientation right of the first index
    double offset = M_PI_2 * octant;
    for (size_t point = 0; point < octSize; ++point)
    {
      double phi = abscissae_[point].phi + offset;
      double theta = abscissae_[point].theta;
      double new_z_value = omegas_[point].z;
      double sinTheta = sqrt(1 - new_z_value * new_z_value);

      new_omega.x = sinTheta * cos(phi); // omegas_[l].x*xsign;
      new_omega.y = sinTheta * sin(phi); // omegas_[l].y*ysign;
      new_omega.z = omegas_[point].z;
      weights_.push_back(weights_[point]);
      omegas_.emplace_back(new_omega);
      abscissae_.emplace_back(phi, theta);
    }
  }

  AngularQuadrature::OptimizeForPolarSymmetry(4.0 * M_PI);
}

void
TriangleQuadrature::FilterMoments(unsigned int scattering_order)
{
  if (m2d_op_built_ and d2m_op_built_ and moments_ >= 0)
  {
    int s_order = static_cast<int>(scattering_order);
    int moments_to_keep = 1 + (s_order * 3 + s_order * s_order) / 2;

    auto m2d_transposed = m2d_op_;
    m_to_ell_em_map_.resize(moments_to_keep);
    std::vector<std::vector<double>> m2dworking;
    std::vector<std::vector<double>> d2mworking;
    for (size_t i{}; i < moments_to_keep; ++i)
    {
      d2mworking.push_back(d2m_op_.at(i));
      m2dworking.push_back(m2d_transposed.at(i));
    }
    m2d_op_ = m2dworking;
    d2m_op_ = d2mworking;
  }
}

void
TriangleQuadrature::MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  int L = static_cast<int>(sn_);
  int L_max = static_cast<int>(scattering_order);
  if (method_ == 0 and m_to_ell_em_map_.empty())
  {
    for (int ell = 0; ell <= L_max; ++ell)
      for (int m = -ell; m <= ell; m += 2)
        m_to_ell_em_map_.emplace_back(ell, m);
  }

  // Standard harmonics
  if (m_to_ell_em_map_.empty() and method_ != 3)
  {
    for (int ell = 0; ell <= L; ++ell)
    {
      for (int m = -ell; m <= ell; m += 2)
      {
        if (ell == L and m >= 0 and ell != 0)
          break;
        else
          m_to_ell_em_map_.emplace_back(ell, m);
      }
    }
  }
  else if (m_to_ell_em_map_.empty())
  {
    for (int ell = 0; ell <= L_max; ++ell)
    {
      for (int m = -ell; m <= ell; m += 2)
      {
        if (ell == L and m >= 0 and ell != 0)
          break;
        else
          m_to_ell_em_map_.emplace_back(ell, m);
      }
    }
  }
}

void
TriangleQuadrature::BuildDiscreteToMomentOperator(unsigned int scattering_order, int dimension)
{
  if (d2m_op_built_)
    return;

  auto inner_product = [](const VecDbl& f, const VecDbl& g, const VecDbl& wt)
  {
    double sum_val = 0.0;
    for (size_t i = 0; i < f.size(); ++i)
      sum_val += f[i] * g[i] * wt[i];
    return sum_val;
  };

  MakeHarmonicIndices(scattering_order, dimension);

  // Standard Sn method
  if (method_ == 0)
  {
    d2m_op_.clear();

    const size_t num_angles = abscissae_.size();
    const size_t num_moms = m_to_ell_em_map_.size();
    for (const auto& ell_em : m_to_ell_em_map_)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);

      for (int n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        cur_mom.push_back(value * weights_[n]);
      }
      d2m_op_.push_back(cur_mom);
    }
    d2m_op_built_ = true;
  }
  else if (method_ == 1)
  {
    d2m_op_.clear();

    if (not m2d_op_built_)
      BuildMomentToDiscreteOperator(scattering_order, dimension);

    const size_t num_angles = abscissae_.size();
    d2m_op_ = Transpose(Inverse(m2d_op_));
    weights_.clear();
    for (const auto& wt : d2m_op_[0])
      weights_.push_back(wt);
    d2m_op_built_ = true;
  }
  if (method_ == 2)
  {
    d2m_op_.clear();

    std::vector<std::vector<double>> cmt;
    unsigned int num_angles = abscissae_.size();
    unsigned int num_moms = m_to_ell_em_map_.size();
    for (const auto& ell_em : m_to_ell_em_map_)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      for (int n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        cur_mom.push_back(value);
      }
      cmt.push_back(cur_mom);
    }

    // Solve for the weights_
    std::vector<double> wt = {4.0 * M_PI};
    for (size_t i = 1; i < m_to_ell_em_map_.size(); ++i)
      wt.emplace_back(0.0);
    weights_ = MatMul(Inverse(cmt), wt);

    for (const auto& ell_em : m_to_ell_em_map_)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      for (int n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        cur_mom.push_back(value * weights_[n]);
      }
      d2m_op_.push_back(cur_mom);
    }
    d2m_op_built_ = true;
  }
  if (method_ == 3) // Using Gram-Schmidt orthogonalization
  {
    d2m_op_.clear();

    // Make the coefficent matrix
    unsigned int num_angles = abscissae_.size();
    unsigned int num_moms = m_to_ell_em_map_.size();
    MatDbl cmt;

    for (const auto& ell_em : m_to_ell_em_map_)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      for (int n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        double value = Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        cur_mom.push_back(value);
      }
      cmt.push_back(cur_mom);
    }

    // Modified Gram-Schmidt
    MatDbl cmt_hat = cmt;
    MatDbl V = cmt_hat;
    MatDbl B = V;
    for (size_t i = 0; i < num_moms; ++i)
    {
      B[i] = VecMul(V[i], 1.0 / inner_product(V[i], V[i], weights_));

      for (size_t k = i + 1; k < num_moms; ++k)
      {
        V[k] = V[k] - VecMul(B[i], inner_product(V[i], V[k], weights_));
      }
    }
    cmt_hat = B;

    // Now normalize the values
    for (int i = 0; i < num_moms; ++i)
    {
      auto ell = m_to_ell_em_map_[i].ell;
      double normal = (4.0 * M_PI) / (2.0 * ell + 1.0);
      double multiplier = sqrt(normal / inner_product(cmt_hat[i], cmt_hat[i], weights_));
      for (int k = 0; k < num_angles; ++k)
        cmt_hat[i][k] *= multiplier;
    }

    // Make the d2m matrix and m2d matrix
    MatDbl holder_m2d;
    for (size_t i = 0; i < num_moms; ++i)
    {
      VecDbl temp_d2m;
      VecDbl temp_m2d;
      auto ell = m_to_ell_em_map_[i].ell;
      for (int k = 0; k < num_angles; ++k)
      {
        temp_m2d.emplace_back(cmt_hat[i][k] * ((2.0 * ell + 1) / (4.0 * M_PI)));
        temp_d2m.emplace_back(cmt_hat[i][k] * weights_[k]);
      }
      d2m_op_.push_back(temp_d2m);
      holder_m2d.push_back(temp_m2d);
    }
    m2d_op_ = holder_m2d;
    d2m_op_built_ = true;
  }

  if (scattering_order < sn_ and method_ != 3 and method_ != 0)
  {
    log.Log0() << "Filtering moements for scattering order " << scattering_order;
    FilterMoments(scattering_order);
  }
}

void
TriangleQuadrature::BuildMomentToDiscreteOperator(unsigned int scattering_order, int dimension)
{
  if (m2d_op_built_)
    return;

  MakeHarmonicIndices(scattering_order, dimension);

  if (method_ == 0)
  {
    const size_t num_angles = abscissae_.size();
    const size_t num_moms = m_to_ell_em_map_.size();

    const auto normalization = std::accumulate(weights_.begin(), weights_.end(), 0.0);

    for (const auto& ell_em : m_to_ell_em_map_)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      for (int n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        double value = ((2.0 * ell_em.ell + 1.0) / normalization) *
                       Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        cur_mom.push_back(value);
      }
      m2d_op_.push_back(cur_mom);
    }
    m2d_op_built_ = true;
  }
  else if (method_ == 1)
  {
    const auto normalization = 4.0 * M_PI;
    const size_t num_angles = abscissae_.size();
    const size_t num_moms = m_to_ell_em_map_.size();
    for (const auto& ell_em : m_to_ell_em_map_)
    {
      double integral = 0.0;
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      for (int n = 0; n < num_angles; ++n)
      {
        const auto& cur_angle = abscissae_[n];
        double value = ((2.0 * ell_em.ell + 1.0) / normalization) *
                       Ylm(ell_em.ell, ell_em.m, cur_angle.phi, cur_angle.theta);
        cur_mom.push_back(value);
      }
      m2d_op_.push_back(cur_mom);
    }
    m2d_op_built_ = true;
  }
  else if (method_ == 2)
  {
    if (not d2m_op_built_)
      BuildDiscreteToMomentOperator(sn_, dimension);
    m2d_op_ = Transpose(Inverse(d2m_op_));
    m2d_op_built_ = true;
  }
  else if (method_ == 3)
  {
    if (not d2m_op_built_)
      BuildDiscreteToMomentOperator(scattering_order, dimension);
    m2d_op_built_ = true;
  }

  if (scattering_order < sn_ and method_ != 3 and method_ != 0)
  {
    log.Log0() << "Filtering moements for scattering order " << scattering_order;
    FilterMoments(scattering_order);
  }
}

} // namespace opensn
