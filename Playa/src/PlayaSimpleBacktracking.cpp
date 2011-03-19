#include "PlayaSimpleBacktracking.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationImpl.hpp"


namespace Playa
{
using std::endl;

SimpleBacktracking::SimpleBacktracking(const ParameterList& params)
  : LineSearchBase(params)
{;}

std::string SimpleBacktracking::description() const
{
  std::ostringstream oss;
  oss << "SimpleBacktracking(maxSteps=" << maxSteps()
      << ", minStepSize=" << minStepSize() << ")";
  return oss.str();
}

LineSearchStatus SimpleBacktracking::search(const RCP<ObjectiveBase>& obj,
  const Vector<double>& x0,
  const double& f0,
  const Vector<double>& direction,
  const double& alphaMax,
  Vector<double>& xn, 
  Vector<double>& gradF,
  double& fVal) const
{
  Tabs tab(0);
  try
  {
    if (verb() > 0)
    {
      Out::root() << tab << "backtracking line search" << endl;
    }
    double alpha = alphaMax;
    if (verb() > 3)
    {
      Out::root() << tab << "line search initial vector " << endl;
      Out::os() << x0 << endl;
      Out::root() << tab << "line search direction " << endl;
      Out::os() << direction << endl;
      Out::root() << tab << "line search initial value " << f0 << endl;
    }

    for (int i = 0; i < maxSteps(); i++)
    {
      Tabs tab1;
      if (verb() > 1)
      {
        Out::root() << tab1 << "line search step " 
                    << i << " alpha=" << alpha << endl;
      }

      xn = x0 + alpha * direction;

      if (verb() > 3)
      {
        Out::root() << tab1 << "line search trial vector " << endl;
        Out::os() << xn << endl;
      }
          
      obj->evalGrad(xn, fVal, gradF);

      if (verb() > 1)
      {
        Out::root() << tab1 << "f=" << fVal << " f0=" << f0 << endl;
      }

      if (alpha < minStepSize()) return LS_StepTooSmall;
          
      if (fVal < f0)
      {
        if (verb() > 0)
        {
          Out::root() << tab1 << "Line search successful: steps = " << i << endl;
        }
        return LS_Success;
      }
      alpha = alpha/2.0;
    }
    return LS_ExceededMaxiters;
  }

  catch(std::exception& e)
  {
    Out::root() << "Exception detected in SimpleBacktracking::search(): " 
                << e.what() << endl;
    return LS_Crashed;
  }
}



}
