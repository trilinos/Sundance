/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnaryFunctor.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceUtils;
using namespace Teuchos;

void UnaryFunctor::eval1(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df_dx) const
{
  evalFDDerivs1(x, nx, f, df_dx);
}


void UnaryFunctor::eval2(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df_dx,
                        double* d2f_dxx) const
{
  evalFDDerivs2(x, nx, f, df_dx, d2f_dxx);
}


void UnaryFunctor::evalFDDerivs1(const double* const x, 
                                 int nx, 
                                 double* f, 
                                 double* df_dx) const
{
  eval0(x, nx, f);

  double w1 = 1.0/h_/12.0;

  for (int i=0; i<nx; i++)
    {
      double xPlus1 = x[i] + h_;
      double xMinus1 = x[i] - h_;
      double fPlus1;
      double fMinus1;
      double xPlus2 = x[i] + 2.0*h_;
      double xMinus2 = x[i] - 2.0*h_;
      double fPlus2;
      double fMinus2;
      eval0(&xPlus1, 1, &fPlus1);
      eval0(&xMinus1, 1, &fMinus1);
      eval0(&xPlus2, 1, &fPlus2);
      eval0(&xMinus2, 1, &fMinus2);
      df_dx[i] = w1*( 8.0*(fPlus1-fMinus1) - (fPlus2 - fMinus2) );
    }
}

void UnaryFunctor::evalFDDerivs2(const double* const x, 
                                 int nx, 
                                 double* f, 
                                 double* df_dx,
                                 double* d2f_dxx) const
{
  eval0(x, nx, f);

  double w1 = 1.0/h_/12.0;
  double w2 = w1/h_;

  for (int i=0; i<nx; i++)
    {
      double xPlus1 = x[i] + h_;
      double xMinus1 = x[i] - h_;
      double fPlus1;
      double fMinus1;
      double xPlus2 = x[i] + 2.0*h_;
      double xMinus2 = x[i] - 2.0*h_;
      double fPlus2;
      double fMinus2;
      eval0(&xPlus1, 1, &fPlus1);
      eval0(&xMinus1, 1, &fMinus1);
      eval0(&xPlus2, 1, &fPlus2);
      eval0(&xMinus2, 1, &fMinus2);
      df_dx[i] = w1*( 8.0*(fPlus1-fMinus1) - (fPlus2 - fMinus2) );
      d2f_dxx[i] = w2*(16.0*(fPlus1 + fMinus1) 
                        - 30.0*f[i] - (fPlus2+fMinus2));
    }
}

bool UnaryFunctor::testDerivs(const double& x, const double& tol) const
{
  Tabs tabs;
  cerr << tabs << endl << tabs 
       << "comparing exact derivs to FD derivs for functor " 
       << name() << endl;

  double fExact;
  double fxExact;
  double fxxExact;

  double fFD;
  double fxFD;
  double fxxFD;

  bool isOK = true;

  /* test first differentiation */
  cerr << tabs << "computing first derivatives at x=" << x << endl;
  {
    Tabs tabs1;

    eval1(&x, 1, &fExact, &fxExact);
    cerr << tabs1 << "Exact: f=" << fExact << " df_dx=" << fxExact << endl; 
    evalFDDerivs1(&x, 1, &fFD, &fxFD);
    cerr << tabs1 << "FD:    f=" << fFD << " df_dx=" << fxFD << endl; 

    double fError = fabs(fFD - fExact)/(fabs(fExact) + h_);
    double fxError = fabs(fxFD - fxExact)/(fabs(fxExact)+h_);
    {
      Tabs tabs2;
      cerr << tabs2 << "| f-f_FD |=" << fError 
           << " | df-df_FD |=" << fxError << endl;
      if (fError > tol) 
        {
          cerr << tabs << "ERROR: function value mismatch!" << endl;
          isOK = false;
        }
      if (fxError > tol) 
        {
          cerr << tabs << "ERROR: first derivative mismatch!" << endl;
          isOK = false;
        }
    }
  }
  
  cerr << tabs << endl;

  
  /* test first and second differentiation */
  cerr << tabs << "computing first and second derivatives at x=" << x << endl;
  {
    Tabs tabs1;

    eval2(&x, 1, &fExact, &fxExact, &fxxExact);
    cerr << tabs1 << "Exact: f=" << fExact 
         << " df_dx=" << fxExact 
         << " d2f_dx2=" << fxxExact 
         << endl; 

    evalFDDerivs2(&x, 1, &fFD, &fxFD, &fxxFD);
    cerr << tabs1 << "FD:    f=" << fFD << " df_dx=" << fxFD  
         << " d2f_dx2=" << fxxFD << endl; 

    double fError = fabs(fFD - fExact)/(fabs(fExact) + h_);
    double fxError = fabs(fxFD - fxExact)/(fabs(fxExact)+h_);
    double fxxError = fabs(fxxFD - fxxExact)/(fabs(fxxExact)+h_);
    {
      Tabs tabs2;
      cerr << tabs2 << "| f-f_FD |=" << fError 
           << " | df-df_FD |=" << fxError 
           << " | d2f-d2f_FD |=" << fxxError << endl;
      if (fError > tol) 
        {
          cerr << tabs << "ERROR: function value mismatch!" << endl;
          isOK = false;
        }
      if (fxError > tol) 
        {
          cerr << tabs << "ERROR: first derivative mismatch!" << endl;
          isOK = false;
        }
      if (fxxError > tol) 
        {
          cerr << tabs << "ERROR: second derivative mismatch!" << endl;
          isOK = false;
        }
    }
  }

  cerr << tabs << endl;
  return isOK;
}

bool UnaryFunctor::testInvalidValue(const double& badValue) const 
{
  Tabs tabs;
  bool detectedError = false;

  cerr << endl << tabs 
       << "testing exception detection for bad value x=" << badValue
       << " for function  " << name() << endl;
  try
    {
      testDerivs(badValue, 1.0);
    }
  catch(std::exception& e)
    {
      detectedError = true;
    }

  if (!detectedError) 
    {
      cerr << tabs << "ERROR: missed detection of an exception at x=" << badValue
           << " for function " << name() << endl;
    }
  else
    {
      cerr << tabs << "the error was detected!" << endl;
    }
  return detectedError;
}


bool UnaryFunctor::test(int nx, const double& tol) const
{
  bool isOK;

  double a = -sqrt(1.5);
  double b = sqrt(2.0);

  /* first test a sample of points in the domain */
  if (domain()->hasLowerBound())
    {
      a = domain()->lowerBound();
    }
  if (domain()->hasUpperBound())
    {
      b = domain()->upperBound();
    }

  double c = (b-a)/((double) nx+1);
  for (int i=1; i<=nx; i++)
    {
      double x = a + i*c;
      if (domain()->hasExcludedPoint() && fabs(x-domain()->excludedPoint())<1.0e-14)
        {
          continue;
        }
      isOK =  testDerivs(x, tol) && isOK ;
    }

  /* Test for exception detection at the excluded point, if any */
  if (domain()->hasExcludedPoint())
    {
      isOK =  testInvalidValue(domain()->excludedPoint()) && isOK ;
    }

  /* Test for exception detection at a point below the lower bound */
  if (domain()->hasLowerBound())
    {
      isOK =  testInvalidValue(domain()->lowerBound() - 0.1) && isOK ;
    }
  
  /* Test for exception detection at a point above the upper bound */
  if (domain()->hasUpperBound())
    {
      isOK =  testInvalidValue(domain()->upperBound() + 0.1) && isOK ;
    }
  
  return isOK;
  

}
