#include "GolubFischer.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <cmath>

/**Master callable function that will return a reference to the abscissae and
 * weights of the discrete angles.*/
AnglePairs& chi_math::GolubFischer::GetDiscreteScatAngles(Tvecdbl& mell)
{
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Getting Discrete Scattering Angles" << '\n';

  for (int m=0; m<mell.size(); m++)
  {
    Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Moment " << m << " " << mell[m];
  }

  std::vector<double> in_mell;
  in_mell = mell;

  int N = in_mell.size()-1;
  int n = (N+1)/2;

  xn_wn_.resize(n, std::pair<double,double>(0.0, 0.0));

  if (N==0)
    return xn_wn_;

  /* Legendre recurrence coefficients */
  Tvecdbl a;   a.resize(2*n,0.0);
  Tvecdbl b;   b.resize(2*n,0.0);
  Tvecdbl c;   c.resize(2*n,0.0);

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "a,b,c:\n";
  for(int j=0; j<2*n; j++)
  {
    a[j] = 0.0;
    b[j] = j/(2.0*j+1);
    c[j] = (j+1.0)/(2*j+1);
    Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2)
      << a[j] << " "
      << b[j] << " "
      << c[j] << " \n";

  }

  MCA(in_mell, a,  b,  c);

  RootsOrtho(n, alpha_, beta_);

  for (int i=0; i<n; i++)
  {
    Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "i " << xn_wn_[i].first << " " << xn_wn_[i].second << '\n';
  }

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Done" << '\n';

  return xn_wn_;
}

/**Applies the Modified Chebyshev Algorithm contained in [1] to find the
 * recursion coefficients for the orthogonal polynomials.*/
void chi_math::GolubFischer::MCA(Tvecdbl& in_mell, Tvecdbl& a, Tvecdbl& b, Tvecdbl& c)
{
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "MCA Start" << '\n';

  int N = in_mell.size()-1;
  int n = (N+1)/2;

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "N " << N << " n " << n << '\n';
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "alpha, beta" << '\n';

  alpha_.resize(n + 1, 0.0);
  beta_.resize(n + 1, 0.0);

  std::vector<std::vector<double>> sigma(n+1,std::vector<double>(2*n+1,0.0));


  for(int ell=0; ell<2*n; ell++)
  {
    sigma[0][ell] = in_mell[ell];
  }

  alpha_[0] = a[0] + c[0] * sigma[0][1] / sigma[0][0];
  beta_[0] = in_mell[0];

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << 0 << " " << alpha_[0] << " " << beta_[0] << "\n";


  for(int k=1; k<n+1; k++)
  {
    for(int ell=k; ell<(2*n-k+1); ell++)
    {
      double sigmakm2ell = 0.0;

      if (k==1)
      {
        sigmakm2ell = 0.0;
      }
      else
      {
        sigmakm2ell = sigma[k-2][ell];
      }
      sigma[k][ell] = c[ell]*sigma[k-1][ell+1]
                      - (alpha_[k - 1] - a[ell]) * sigma[k - 1][ell]
                      - beta_[k - 1] * sigmakm2ell
                      +b[ell]*sigma[k-1][ell-1];
    }
    alpha_[k] = a[k]
                -c[k-1]*(sigma[k-1][k]/sigma[k-1][k-1])
                +c[k]*(sigma[k][k+1]/sigma[k][k]);
    beta_[k] = c[k - 1] * sigma[k][k] / sigma[k - 1][k - 1];

    Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << k << " " << alpha_[k] << " " << beta_[k] << "\n";
  }

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Done" << '\n';
}

/**Finds the roots of the orthogonal polynomial.*/
void chi_math::GolubFischer::RootsOrtho(int& N, Tvecdbl& in_alpha, Tvecdbl& in_beta)
{
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "RootsOrtho Start" << '\n';

  int maxiters = 1000;
  double tol = 1.0e-6;
  double adder = 0.999*2/std::max(N-1,1);

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Check 1: Init guess" << '\n';

  Tvecdbl xn; xn.resize(N, 0.0);
  Tvecdbl wn; wn.resize(N, 0.0);

  for(int i=0; i<N; i++)
  {
    xn[i] = -0.999+i*adder;
    Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "x[" << i << "]=" << xn[i] << "\n";
  }



  Tvecdbl norm; norm.resize(N+1, 0.0);
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Check 2 " << in_beta[0] << '\n';
  norm[0] = in_beta[0];
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Check 3a norms" << '\n';
  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << norm[0] << '\n';
  for (int i=1; i<(N+1); i++)
  {
    norm[i]=in_beta[i]*norm[i-1];
    Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << norm[i] << '\n';
  }

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Check 3" << '\n';

  for (int k=0; k<N; k++)
  {
    int i = 0;

    while (i<maxiters)
    {
      double xold = xn[k];
      double a = Ortho(N, xold, in_alpha, in_beta);
      double b = dOrtho(N, xold, in_alpha, in_beta);
      double c = 0;

      for (int j=0; j<k; j++)
      {
        c = c+(1.0/(xold-xn[j]));
      }

      double xnew = xold-(a/(b-a*c));
      if (std::isnan(xnew))
      {
        Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "xnew " << i << " " << xnew << " y=" << a << std::endl;
        Chi::Exit(EXIT_FAILURE);
      }

      double res = std::fabs(xnew-xold);
      xn[k] = xnew;

      Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "xnew " << i << " " << xnew << " y=" << a << std::endl;

      if (res<tol)
      {
        break;
      }

      i = i+1;
    }//while

  }//for k

  for (int i=0; i<N-1; i++)
  {
    for (int j=0; j<N-i-1; j++)
    {
      if (xn[j]>xn[j+1])
      {
        double tempx = xn[j+1];
        double tempw = wn[j+1];
        xn[j+1] = xn[j];
        wn[j+1] = wn[j];
        xn[j] = tempx;
        wn[j] = tempw;
      }
    }
  }

  for (int i=0; i<N; i++)
  {
    wn[i] = 0.0;

    for (int k=0; k<N; k++)
    {
      wn[i] += Ortho(k, xn[i], in_alpha, in_beta)*
               Ortho(k, xn[i], in_alpha, in_beta)/norm[k];
    }

    wn[i] = 1.0/wn[i];
  }
  for (int i=0; i<N; i++)
  {
    xn_wn_[i].first = xn[i];
    xn_wn_[i].second = wn[i];
  }

  Chi::log.Log(chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2) << "Done" << '\n';
}

/**Computes the function evaluation of the orthogonal polynomials.*/
double chi_math::GolubFischer::Ortho(int ell, double x, Tvecdbl& in_alpha, Tvecdbl& in_beta)
{
  if (ell==0)
  {
    return 1;
  }

  if (ell==1)
  {
    return (x-in_alpha[0]);
  }

  double Pnm1 = 1.0;
  double Pn   = (x-in_alpha[0]);
  double Pnp1 = 0.0;

  for (int n=2; n<ell+1; n++)
  {
    int ns = n-1;
    Pnp1 = (x-in_alpha[ns])*Pn-in_beta[ns]*Pnm1;
    Pnm1 = Pn;
    Pn = Pnp1;
  }

  return Pnp1;
}

/**Computes the derivative of the orthogonal polynomials.*/
double chi_math::GolubFischer::dOrtho(int ell, double x, Tvecdbl& in_alpha, Tvecdbl& in_beta)
{

  double eps = 0.000001;
  double y2 = Ortho(ell, x+eps, in_alpha, in_beta);
  double y1 = Ortho(ell, x-eps, in_alpha, in_beta);

  double m = (y2-y1)/2.0/eps;

  return m;
}