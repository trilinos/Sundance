/* @HEADER@ */
/* @HEADER@ */

#include "SundanceStdMathOps.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


PowerFunctor::PowerFunctor(const double& p) 
  : UnaryFunctor("pow("+Teuchos::toString(p)+")"), p_(p)
{;}

void PowerFunctor::eval(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df) const
{
  if (fabs(p_) > 1.0e-16)
    {
      for (int i=0; i<nx; i++) 
        {
          double px = ::pow(x[i], p_-1);
          df[i] = p_*px;
          f[i] = x[i]*px;
        }
    }
  else
    {
      for (int i=0; i<nx; i++) 
        {
          f[i] = 1.0;
          df[i] = 0.0;
        }
    }
}

void PowerFunctor::eval(const double* const x, 
                        int nx, 
                        double* f) const
{
  for (int i=0; i<nx; i++) 
    {
      f[i] = ::pow(x[i], p_);
    }
}
