/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaRosenbrock.hpp"

namespace Playa
{

void Rosenbrock::eval(const Vector<double>& x, double& f) const
{
  f = 0.0;
  for (int i=0; i<N_; i++)
  {
    double p = x[2*i+1] - x[2*i]*x[2*i];
    double q = 1.0-x[2*i];
    f += alpha_ * p*p + q*q;
  }
}

void Rosenbrock::evalGrad(const Vector<double>& x, double& f,
  Vector<double>& df) const
{
  f = 0.0;
  df.zero();
  for (int i=0; i<N_; i++)
  {
    double p = x[2*i+1] - x[2*i]*x[2*i];
    double q = 1.0-x[2*i];
    f += alpha_ * p*p + q*q;
    df[2*i] += -4.0*alpha_*p*x[2*i] - 2.0*q;
    df[2*i+1] += 2.0*alpha_*p; 
  }  
}

Vector<double> Rosenbrock::getInit() const
{
  Vector<double> rtn = vs_.createMember();
  rtn.setToConstant(-1.0);
  return rtn;
}


}
