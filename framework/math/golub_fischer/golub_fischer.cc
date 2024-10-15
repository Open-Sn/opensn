// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/golub_fischer/golub_fischer.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <cmath>

namespace opensn
{

std::vector<std::pair<double, double>>&
GolubFischer::GetDiscreteScatAngles(std::vector<double>& mell)
{
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Getting Discrete Scattering Angles" << '\n';

  for (int m = 0; m < mell.size(); ++m)
  {
    log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Moment " << m << " " << mell[m];
  }

  std::vector<double> in_mell;
  in_mell = mell;

  int N = in_mell.size() - 1;
  int n = (N + 1) / 2;

  xn_wn_.resize(n, std::pair<double, double>(0.0, 0.0));

  if (N == 0)
    return xn_wn_;

  /* Legendre recurrence coefficients */
  std::vector<double> a;
  a.resize(2 * n, 0.0);
  std::vector<double> b;
  b.resize(2 * n, 0.0);
  std::vector<double> c;
  c.resize(2 * n, 0.0);

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "a,b,c:\n";
  for (int j = 0; j < 2 * n; ++j)
  {
    a[j] = 0.0;
    b[j] = j / (2.0 * j + 1);
    c[j] = (j + 1.0) / (2 * j + 1);
    log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << a[j] << " " << b[j] << " " << c[j] << " \n";
  }

  MCA(in_mell, a, b, c);

  RootsOrtho(n, alpha_, beta_);

  for (int i = 0; i < n; ++i)
  {
    log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2)
      << "i " << xn_wn_[i].first << " " << xn_wn_[i].second << '\n';
  }

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Done" << '\n';

  return xn_wn_;
}

void
GolubFischer::MCA(std::vector<double>& mell,
                  std::vector<double>& a,
                  std::vector<double>& b,
                  std::vector<double>& c)
{
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "MCA Start" << '\n';

  int N = mell.size() - 1;
  int n = (N + 1) / 2;

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "N " << N << " n " << n << '\n';
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "alpha, beta" << '\n';

  alpha_.resize(n + 1, 0.0);
  beta_.resize(n + 1, 0.0);

  std::vector<std::vector<double>> sigma(n + 1, std::vector<double>(2 * n + 1, 0.0));

  for (int ell = 0; ell < 2 * n; ++ell)
  {
    sigma[0][ell] = mell[ell];
  }

  alpha_[0] = a[0] + c[0] * sigma[0][1] / sigma[0][0];
  beta_[0] = mell[0];

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << 0 << " " << alpha_[0] << " " << beta_[0] << "\n";

  for (int k = 1; k < n + 1; ++k)
  {
    for (int ell = k; ell < (2 * n - k + 1); ++ell)
    {
      double sigmakm2ell = 0.0;

      if (k == 1)
      {
        sigmakm2ell = 0.0;
      }
      else
      {
        sigmakm2ell = sigma[k - 2][ell];
      }
      sigma[k][ell] = c[ell] * sigma[k - 1][ell + 1] -
                      (alpha_[k - 1] - a[ell]) * sigma[k - 1][ell] - beta_[k - 1] * sigmakm2ell +
                      b[ell] * sigma[k - 1][ell - 1];
    }
    alpha_[k] = a[k] - c[k - 1] * (sigma[k - 1][k] / sigma[k - 1][k - 1]) +
                c[k] * (sigma[k][k + 1] / sigma[k][k]);
    beta_[k] = c[k - 1] * sigma[k][k] / sigma[k - 1][k - 1];

    log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << k << " " << alpha_[k] << " " << beta_[k] << "\n";
  }

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Done" << '\n';
}

void
GolubFischer::RootsOrtho(int& N, std::vector<double>& alpha, std::vector<double>& beta)
{
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "RootsOrtho Start" << '\n';

  int maxiters = 1000;
  double tol = 1.0e-6;
  double adder = 0.999 * 2 / std::max(N - 1, 1);

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Check 1: Init guess" << '\n';

  std::vector<double> xn;
  xn.resize(N, 0.0);
  std::vector<double> wn;
  wn.resize(N, 0.0);

  for (int i = 0; i < N; ++i)
  {
    xn[i] = -0.999 + i * adder;
    log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "x[" << i << "]=" << xn[i] << "\n";
  }

  std::vector<double> norm;
  norm.resize(N + 1, 0.0);
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Check 2 " << beta[0] << '\n';
  norm[0] = beta[0];
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Check 3a norms" << '\n';
  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << norm[0] << '\n';
  for (int i = 1; i < (N + 1); ++i)
  {
    norm[i] = beta[i] * norm[i - 1];
    log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << norm[i] << '\n';
  }

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Check 3" << '\n';

  for (int k = 0; k < N; ++k)
  {
    int i = 0;

    while (i < maxiters)
    {
      double xold = xn[k];
      double a = Ortho(N, xold, alpha, beta);
      double b = dOrtho(N, xold, alpha, beta);
      double c = 0;

      for (int j = 0; j < k; ++j)
      {
        c = c + (1.0 / (xold - xn[j]));
      }

      double xnew = xold - (a / (b - a * c));
      if (std::isnan(xnew))
      {
        log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2)
          << "xnew " << i << " " << xnew << " y=" << a << std::endl;
        Exit(EXIT_FAILURE);
      }

      double res = std::fabs(xnew - xold);
      xn[k] = xnew;

      log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2)
        << "xnew " << i << " " << xnew << " y=" << a << std::endl;

      if (res < tol)
      {
        break;
      }

      i = i + 1;
    } // while

  } // for k

  for (int i = 0; i < N - 1; ++i)
  {
    for (int j = 0; j < N - i - 1; ++j)
    {
      if (xn[j] > xn[j + 1])
      {
        double tempx = xn[j + 1];
        double tempw = wn[j + 1];
        xn[j + 1] = xn[j];
        wn[j + 1] = wn[j];
        xn[j] = tempx;
        wn[j] = tempw;
      }
    }
  }

  for (int i = 0; i < N; ++i)
  {
    wn[i] = 0.0;

    for (int k = 0; k < N; ++k)
    {
      wn[i] += Ortho(k, xn[i], alpha, beta) * Ortho(k, xn[i], alpha, beta) / norm[k];
    }

    wn[i] = 1.0 / wn[i];
  }
  for (int i = 0; i < N; ++i)
  {
    xn_wn_[i].first = xn[i];
    xn_wn_[i].second = wn[i];
  }

  log.Log(Logger::LOG_LVL::LOG_0VERBOSE_2) << "Done" << '\n';
}

double
GolubFischer::Ortho(int ell, double x, std::vector<double>& alpha, std::vector<double>& beta)
{
  if (ell == 0)
  {
    return 1;
  }

  if (ell == 1)
  {
    return (x - alpha[0]);
  }

  double Pnm1 = 1.0;
  double Pn = (x - alpha[0]);
  double Pnp1 = 0.0;

  for (int n = 2; n < ell + 1; ++n)
  {
    int ns = n - 1;
    Pnp1 = (x - alpha[ns]) * Pn - beta[ns] * Pnm1;
    Pnm1 = Pn;
    Pn = Pnp1;
  }

  return Pnp1;
}

double
GolubFischer::dOrtho(int ell, double x, std::vector<double>& alpha, std::vector<double>& beta)
{

  double eps = 0.000001;
  double y2 = Ortho(ell, x + eps, alpha, beta);
  double y1 = Ortho(ell, x - eps, alpha, beta);

  double m = (y2 - y1) / 2.0 / eps;

  return m;
}

} // namespace opensn
