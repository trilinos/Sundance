#include "PlayaBasicLMBFGS.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLineSearchBuilder.hpp"
#include "PlayaOptConvergenceTestBuilder.hpp"

namespace Playa
{
using std::endl;

BasicLMBFGS::BasicLMBFGS(
  const ParameterList& params
  )
  : LineSearchBasedOptBase(params),
    memSize_(getParameter<int>(params, "Max Memory Size"))
{}


RCP<DirectionGeneratorBase> 
BasicLMBFGS::makeDirectionGenerator() const 
{
  return rcp(new BasicLMBFGSDirection(memSize_));
}


BasicLMBFGSDirection::BasicLMBFGSDirection(int memSize)
  : memSize_(memSize),
    xPrev_(),    
    gradPrev_(),
    sMem_(),
    yMem_()
{}

bool BasicLMBFGSDirection::generateDirection(
  const RCP<ObjectiveBase>& obj,
  const Vector<double>& xCur,
  const Vector<double>& gradCur,
  const double& fCur,
  Vector<double>& p)
{
  Vector<double> q = gradCur.copy();
  int numStored = sMem_.size();
  Array<double> alpha(numStored);
  Array<double> rho(numStored);

  for (int i=numStored-1; i>=0; i--)
  {
    rho[i] = 1.0/(sMem_[i]*yMem_[i]);
    alpha[i] = rho[i] * (sMem_[i]*q);
    q = q - alpha[i]*yMem_[i];
  }
  
  double gamma;
  if (numStored > 0)
  {
    int j = numStored-1;
    gamma = (sMem_[j]*yMem_[j])/(yMem_[j]*yMem_[j]);
  }
  else
  {
    gamma = obj->getInvHScale();
  }

  Vector<double> r = gamma*q;

  for (int i=0; i<numStored; i++)
  {
    double beta = rho[i]*(yMem_[i]*r);
    r = r + (alpha[i]-beta)*sMem_[i];
  }

  p = -1.0*r;

  if (xPrev_.ptr().get() != 0)
  {
    Vector<double> s = xCur - xPrev_;
    Vector<double> y = gradCur - gradPrev_;
    sMem_.push_back(s);
    yMem_.push_back(y);
    if ((int) sMem_.size() > memSize_)
    {
      sMem_.pop_front();
      yMem_.pop_front();
    }
  }
  
  xPrev_.acceptCopyOf(xCur);
  gradPrev_.acceptCopyOf(gradCur);

  return true;
}


}
